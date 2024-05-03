using Flux, CUDA, Statistics, ProgressMeter, Plots, BSON

# Brooks Corey relperm formula that we are trying to learn

# s = saturation (variable range 0.0, 1.0)
# n = Exponents for each phase (range from 1.0 to 6.0) set to constant 2
# sr = Residual saturations for each phase (set to 0.2)
# kwm = The maximum relative permeabilities (range from 0.0 to 1.0) set to constant 1.0
# sr_tot = Total residual saturation over all phases i.e. S_or + S_gr + S_wr. (range 0.0 to 1.0) should sr*number of phases

# the output should be between 0 and 1

function brooks_corey_relperm(s::T, n::Real, sr::Real, kwm::Real, sr_tot::Real) where T
    den = 1 - sr_tot
    sat = (s - sr) / den
    sat = clamp(sat, zero(T), one(T))
    return kwm*sat^n
end

# Generate some data for training a model to represent BrooksCorey function
#training_sat = rand(Float32, 1000);                                    # 1×1000 Matrix{Float32} [0.0,1.0]
training_sat = collect(range(Float32(0), stop=Float32(1), length=10000))
rel_perm_analytical = Array{Float32, 1}(undef, 10000); # Creates a one-dimensional array of Float32 with 1000 elements

training_sat = reshape(training_sat, 1, :)
rel_perm_analytical = reshape(rel_perm_analytical, 1, :)

println("Size of training_sat: ", size(training_sat))
println("Size of rel_perm_analytical: ", size(rel_perm_analytical))
println("Shape of input data: ", size(training_sat))
println("Shape of output data: ", size(rel_perm_analytical))
println("Number of training samples: ", length(training_sat))

for i in eachindex(training_sat)
    rel_perm_analytical[i] =  brooks_corey_relperm(training_sat[i], 2.0, 0.2, 1.0, 0.4)
end

plot(vec(training_sat), vec(rel_perm_analytical), label="Brooks-Corey RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")

# Define our model, a multi-layer perceptron with two hidden layers of size 100
model = Chain(
    Dense(1 => 100, relu; init=Flux.glorot_normal),   # activation function inside layer
    Dense(100 => 100, relu; init=Flux.glorot_normal),
    Dense(100 => 100, relu; init=Flux.glorot_normal),
    Dense(100 => 1; init=Flux.glorot_normal),
    relu) |> gpu        # move model to GPU, if available

# The model takes in the saturation with the shape (1xN)
rel_perm_predicted = model(training_sat |> gpu) |> cpu                                 # 1×1000 Matrix{Float32}
# The output of hte model is the relative Permeability with shape (1xN)

# To train the model, we use batches of 64 samples
loader = Flux.DataLoader((training_sat, rel_perm_analytical) |> gpu, batchsize=256, shuffle=true);

optim = Flux.setup(Flux.Adam(0.0001), model)  # will store optimiser momentum, etc.

# Training loop, using the whole data set 1000 times:
losses = []
@showprogress for epoch in 1:1000
    for (x, y) in loader
        loss, grads = Flux.withgradient(model) do m
            # Evaluate model and loss inside gradient context:
            y_hat = m(x)
            Flux.mse(y_hat, y)
        end
        Flux.update!(optim, model, grads[1])
        push!(losses, loss)  # logging, outside gradient context
    end
end

# plot loss function

plot(losses; xaxis=(:log10, "iteration"),
    yaxis="loss", label="per batch")
n = length(loader)
plot!(n:n:length(losses), mean.(Iterators.partition(losses, n)),
    label="epoch mean", dpi=200)

# Predict on the trained model
rel_perm_pred = model(training_sat |> gpu) |> cpu

plot(vec(training_sat), vec(rel_perm_analytical), label="Brooks-Corey RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")
plot!(vec(training_sat), vec(rel_perm_pred), label="ML model RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")

using BSON: @save

@save "BrooksCoreyMLModel.bson" model

using Flux, BSON

BSON.@load "BrooksCoreyMLModel.bson" model

# test on random inputs, different to the training set

# Generate 1000 random numbers between 0 and 1
testing_sat = rand(Float32, 1000)
# sort for easier plotting
testing_sat = sort(testing_sat)
testing_sat = reshape(testing_sat, 1, :)

# Calculate analytical solution using Brooks Corey relperm
test_y = Array{Float32, 1}(undef, 1000)
test_y = reshape(test_y, 1, :)
for i in eachindex(testing_sat)
    test_y[i] = brooks_corey_relperm(testing_sat[i], 2.0, 0.2, 1.0, 0.4)
end

# Predict on the trained model
pred_y = model(testing_sat |> gpu) |> cpu

plot(vec(testing_sat), vec(test_y), label="Brooks-Corey RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")
plot!(vec(testing_sat), vec(pred_y), label="ML model RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")