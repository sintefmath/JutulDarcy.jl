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
#training_sat = rand(Float64, 1000);                                    # 1×1000 Matrix{Float64} [0.0,1.0]
training_sat = collect(range(Float64(0), stop=Float64(1), length=500000))
#training_sat = vcat(zeros(1000), training_sat, ones(1000)) # important to be well behaved at 0.0 and 1.0. so adding training data in this range
rel_perm_analytical = Array{Float64, 1}(undef, 500000);

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

struct CustomModel{T <: Chain} # Parameter to avoid type instability
    chain::T
  end

function (m::CustomModel)(x)
    return (1.0.-x).*0.5.*(tanh.(m.chain(x)) .+ 1.0).*x .+ x
    #min.(max.(0.0, m.chain(x)), 1.0) # clamping between 0 and 1
  end

# Call @layer to allow for training. Described below in more detail.
Flux.@layer CustomModel


# The model takes in the saturation with the shape (1xN)
# The output of the model is the relative Permeability with shape (1xN)


# Define our model, a multi-layer perceptron with three hidden layers of size 50
# use tanh as we want a smooth first derivative
MLP = f64(Chain( # use float64 to match Jutul
    Dense(1 => 50, tanh; init=Flux.glorot_normal),   # activation function inside layer
    Dense(50 => 50, tanh; init=Flux.glorot_normal),
    Dense(50 => 50, tanh; init=Flux.glorot_normal),
    Dense(50 => 50, tanh; init=Flux.glorot_normal),
    # use sigmoid in final activation to ensure output is between 0 and 1
    Dense(50 => 1, sigmoid; init=Flux.glorot_normal)))

    BrooksCoreyMLModel = f64(MLP |> gpu) # move model to GPU, if available)

    # alternatively use custom final activation function
    #BrooksCoreyMLModel = CustomModel(MLP)


# To train the model, we use batches of 10000
loader = Flux.DataLoader((training_sat, rel_perm_analytical) |> gpu, batchsize=10000, shuffle=true);

optim = Flux.setup(Flux.Adam(0.0001), BrooksCoreyMLModel)  # will store optimiser momentum, etc.

# Training loop, using the whole data set 1000 times:
losses = []
@showprogress for epoch in 1:1000
    for (x, y) in loader
        loss, grads = Flux.withgradient(BrooksCoreyMLModel) do m
            # Evaluate model and loss inside gradient context:
            y_hat = m(x)
            Flux.mse(y_hat, y)
        end
        Flux.update!(optim, BrooksCoreyMLModel, grads[1])
        push!(losses, loss)  # logging, outside gradient context
    end
end


# plot loss function

plot(losses; xaxis=(:log10, "iteration"),
    yaxis=(:log10, "loss"), label="per batch")
n = length(loader)
plot!(n:n:length(losses), mean.(Iterators.partition(losses, n)),
    label="epoch mean", dpi=200)

# Predict on the trained model
rel_perm_pred = BrooksCoreyMLModel(training_sat |> gpu) |> cpu

plot(vec(training_sat), vec(rel_perm_analytical), label="Brooks-Corey RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")
plot!(vec(training_sat), vec(rel_perm_pred), label="ML model RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")

plot(vec(training_sat), vec(rel_perm_pred - rel_perm_analytical), label="ML model RelPerm - analytical", xlabel="Saturation", ylabel="Relative Permeability", title="Relative Permeability Error")

println("Norm: ", norm(rel_perm_pred - rel_perm_analytical))

BrooksCoreyMLModel = f64(BrooksCoreyMLModel |> cpu)
@save "BrooksCoreyMLModel.bson" BrooksCoreyMLModel


BSON.@load "BrooksCoreyMLModel.bson" BrooksCoreyMLModel
BrooksCoreyMLModel = f64(BrooksCoreyMLModel |> gpu)

# test on random inputs, different to the training set

# Generate 1000 random numbers between 0 and 1
testing_sat = rand(Float64, 1000)
# add 0 and 1 to testing_sat
testing_sat[1] = 0.0
testing_sat[end] = 1.0
# sort for easier plotting
testing_sat = sort(testing_sat)
testing_sat = reshape(testing_sat, 1, :)

# Calculate analytical solution using Brooks Corey relperm
test_y = Array{Float64, 1}(undef, 1000)
test_y = reshape(test_y, 1, :)
for i in eachindex(testing_sat)
    test_y[i] = brooks_corey_relperm(testing_sat[i], 2.0, 0.2, 1.0, 0.4)
end

# Predict on the trained model
pred_y = BrooksCoreyMLModel(testing_sat |> gpu) |> cpu

plot(vec(testing_sat), vec(test_y), label="Brooks-Corey RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")
plot!(vec(testing_sat), vec(pred_y), label="ML model RelPerm", xlabel="Saturation", ylabel="Relative Permeability", title="Saturation vs. Relative Permeability")

function f_brooks_corey_relperm(s::T, n::Real=2, sr::Real=0.2, kwm::Real=1.0, sr_tot::Real=0.4) where T
    den = 1 - sr_tot
    sat = (s - sr) / den
    sat = clamp(sat, zero(T), one(T))
    return kwm*sat^n
end

#=
Evaluate the gradient of the Brooks Corey relperm function and ML model


For the parameters we have used, we can easily calculate the analytical expression of the gradient for our function (ignoring the clamp function for now):

function f_brooks_corey_relperm(s)
den = 0.6
sat = kwm*(s - sr) / den = (1*(s - 0.2) / 0.6)
return 1*sat^2

->

f(s) = 1*(s/0.6 - 0.2/0.6)^2

f'(s) = 2*((s - 0.2)/0.6)*1/0.6  = 50/9 * (s - 0.2)
=#

using ForwardDiff

BrooksCoreyMLModel = BrooksCoreyMLModel |> cpu

println("f_brooks_corey_relperm(0.5): ", f_brooks_corey_relperm(0.5))

s = 0.5

ForwardDiff.derivative(f_brooks_corey_relperm, s)

f_gradients =  Array{Float64, 1}(undef, 1000);

for i in eachindex(testing_sat)
    f_gradients[i] = ForwardDiff.derivative(f_brooks_corey_relperm, testing_sat[i])
end

f_model_gradients_zygote = gradient(testing_sat -> sum(BrooksCoreyMLModel(testing_sat)), testing_sat)

f_model_gradients_forwardDiff = ForwardDiff.gradient(testing_sat -> sum(BrooksCoreyMLModel(testing_sat)), testing_sat)
#f_model_gradients_forwardDiff = diag(ForwardDiff.jacobian(testing_sat -> BrooksCoreyMLModel(testing_sat), testing_sat))

f_gradients_analytic =  Array{Float64, 1}(undef, 1000);
for i in eachindex(testing_sat)
    f_gradients_analytic[i] = 50/9*(testing_sat[i]-0.2)
    # manually adding the effect of the clamp function
    if (testing_sat[i] < 0.2) || (testing_sat[i] > 0.8)
        f_gradients_analytic[i] = 0
    end
end

println("Every 100th element of f_gradients_analytic:")
println(f_gradients_analytic[1:100:end])

plot(vec(testing_sat[1:1000]), vec(f_gradients[1:1000]), marker=(:circle,2), label="Brooks Coorey function derivative", xlabel="Saturation", ylabel="Relative Permeability derivative", title="Derivative comparison")
plot!(vec(testing_sat[1:1000]), vec(f_gradients_analytic[1:1000]), marker=(:circle,2), label="Brooks Coorey analytical derivative", xlabel="Saturation", ylabel="Relative Permeability", title="Derivative comparison")
plot!(vec(testing_sat[1:1000]), vec(f_model_gradients_zygote[1][1:1000]), marker=(:circle,2), label="Brooks Coorey ML model derivative (Zygote)", xlabel="Saturation", ylabel="Relative Permeability", title="Derivative comparison")
plot!(vec(testing_sat[1:1000]), vec(f_model_gradients_forwardDiff[1:1000]), marker=(:circle,2), label="Brooks Coorey ML model derivative (ForwardDiff)", xlabel="Saturation", ylabel="Relative Permeability", title="Derivative comparison")

plot(vec(testing_sat), vec(vec(f_model_gradients_forwardDiff) - f_gradients), label="ML model derivative - analytical", xlabel="Saturation", ylabel="Relative Permeability derivative error", title="Derivative error")

plot(vec(testing_sat), vec(vec(f_model_gradients_zygote[1]) - vec(f_model_gradients_forwardDiff)), label="Zygote - ForwardDiff", xlabel="Saturation", ylabel="Relative Permeability derivative", title="Zygote - ForwardDiff")

println("f_brooks_corey_relperm(0.0):")
print(f_brooks_corey_relperm(0.0))
s = 0.0
println("\nDerivative at s = 0.0:")
println(ForwardDiff.derivative(f_brooks_corey_relperm, s))

println("\nf_brooks_corey_relperm(1.0):")
print(f_brooks_corey_relperm(1.0))
s = 1.0
println("\nDerivative at s = 1.0:")
println(ForwardDiff.derivative(f_brooks_corey_relperm, s))

println("\nf_brooks_corey_relperm(0.5):")
print(f_brooks_corey_relperm(0.5))
s = 0.5
println("\nDerivative at s = 0.5:")
println(ForwardDiff.derivative(f_brooks_corey_relperm, s))

s = [0.5]
pred_y_0 = BrooksCoreyMLModel(s)
println("\nBrooksCoreyMLModel([0.5]):")
println(pred_y_0)

s = [0.5]
println("\nGradient of BrooksCoreyMLModel at s = [0.5]:")
println(ForwardDiff.gradient(s -> BrooksCoreyMLModel(s)[1], s))

s = [0.0]
pred_y_0 = BrooksCoreyMLModel(s)
println("\nBrooksCoreyMLModel([0.0]):")
println(pred_y_0)

s = [0.0]
println("\nGradient of BrooksCoreyMLModel at s = [0.0]:")
println(ForwardDiff.gradient(s -> BrooksCoreyMLModel(s)[1], s))

s = [1.0]
pred_y_0 = BrooksCoreyMLModel(s)
println("\nBrooksCoreyMLModel([1.0]):")
println(pred_y_0)

s = [1.0]
println("\nGradient of BrooksCoreyMLModel at s = [1.0]:")
println(ForwardDiff.gradient(s -> BrooksCoreyMLModel(s)[1], s))