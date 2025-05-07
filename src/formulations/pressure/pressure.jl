struct PressureFormulation <: Jutul.JutulFormulation

end

const PressureModel = SimulationModel{<:Any, <:Any, PressureFormulation, <:Any}

struct PressureReductionFactors <: Jutul.VectorVariables end

Jutul.values_per_entity(model::PressureModel, ::PressureReductionFactors) = number_of_components(model.system)

struct PressureEquation{T} <: JutulEquation where T<:ConservationLaw
    conservation::T
end

struct PressureEquationTPFAStorage{A, HC, B}
    accumulation::A
    accumulation_symbol::Symbol
    half_face_flux_cells::HC
    buf::B
end

function PressureEquationTPFAStorage(model, eq::PressureEquation; kwarg...)
    ceq = eq.conservation
    ceq.flow_discretization::TwoPointPotentialFlowHardCoded
    number_of_equations = 1
    D, ctx = model.domain, model.context
    cell_entity = Cells()
    face_entity = Faces()
    nc = count_active_entities(D, cell_entity, for_variables = false)
    nf = count_active_entities(D, face_entity, for_variables = false)
    nhf = number_of_half_faces(ceq.flow_discretization)
    face_partials = degrees_of_freedom_per_entity(model, face_entity)
    @assert face_partials == 0 "Only supported for cell-centered discretization"
    alloc = (n, entity, n_entities_pos) -> CompactAutoDiffCache(number_of_equations, n, model,
                                                                                entity = entity, n_entities_pos = n_entities_pos, 
                                                                                context = ctx; kwarg...)
    # Accumulation terms
    acc = alloc(nc, cell_entity, nc)
    # Source terms - as sparse matrix
    t_acc = eltype(acc.entries)
    # src = sparse(zeros(0), zeros(0), zeros(t_acc, 0), size(acc.entries)...)
    # Half face fluxes - differentiated with respect to pairs of cells
    hf_cells = alloc(nhf, cell_entity, nhf)
    # # Half face fluxes - differentiated with respect to the faces
    # if face_partials > 0
    #     hf_faces = alloc(nf, face_entity, nhf)
    # else
    #     hf_faces = nothing
    # end
    buf = zeros(t_acc, number_of_components(model.system))
    return PressureEquationTPFAStorage(acc, conserved_symbol(ceq), hf_cells, buf)
end

function Jutul.setup_equation_storage(model,
        eq::PressureEquation{ConservationLaw{A, B, C, D}}, storage; extra_sparsity = nothing, kwarg...
        ) where {A, B<:TwoPointPotentialFlowHardCoded, C, D}
    return PressureEquationTPFAStorage(model, eq; kwarg...)
end

Jutul.discretization(eq::PressureEquation) = Jutul.discretization(eq.conservation)

include("variables.jl")
include("overloads.jl")
include("functions.jl")
include("bc_and_sources.jl")
include("multimodel/multimodel_pressure.jl")