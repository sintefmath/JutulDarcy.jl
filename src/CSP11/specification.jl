# SPE11 statement: Tables 1, 2, 4 and 5; Sections 2.2–2.6, 3.1–3.7 and 4.2–4.7.
const spe11_year = 365*si_unit(:day)
const spe11_gravity = 9.81
const spe11_size_x = 8400.0
const spe11_size_y = 5000.0
const spe11_size_z = 1200.0
const spe11_geometry_scale = (3000.0, 1000.0)
const spe11_buffer_length = 5e4
# Storage regularization for facies 7. The statement specifies zero porosity.
const spe11_facies7_porosity = 0.001
const spe11_wetting_immobile_saturation = (0.32, 0.14, 0.12, 0.12, 0.12, 0.10)
const spe11_nonwetting_immobile_saturation = 0.10
const spe11_leverett_coefficient = 6.12e-3

function check_spe11_case(case)
    case in (:a, :b, :c) || throw(ArgumentError("Only SPE11 cases :a, :b and :c are supported"))
    return case
end

function spe11_case_specification(case)
    check_spe11_case(case)
    laboratory = case == :a
    three_dimensional = case == :c
    permeability = laboratory ? (4e-11, 5e-10, 1e-9, 2e-9, 4e-9, 1e-8, 0.0) :
        (1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12, 0.0)
    porosity = laboratory ? (0.44, 0.43, 0.44, 0.45, 0.43, 0.46, spe11_facies7_porosity) :
        (0.10, 0.20, 0.20, 0.20, 0.25, 0.35, spe11_facies7_porosity)
    entry_pressure = laboratory ? (1500.0, 300.0, 100.0, 25.0, 10.0, 1.0) :
        ntuple(i -> spe11_leverett_coefficient*sqrt(porosity[i]/permeability[i]), 6)
    boxes = laboratory ? ((1.1, 2.8, 0.0, 0.6), (0.0, 1.1, 0.6, 1.2), (1.1, 2.6, 0.1, 0.4)) :
        three_dimensional ? ((3300.0, 8300.0, 0.0, 750.0), (100.0, 3300.0, 750.0, 1350.0), (3300.0, 7800.0, 250.0, 550.0)) :
        ((3300.0, 8300.0, 0.0, 600.0), (100.0, 3300.0, 600.0, 1200.0), (3300.0, 7800.0, 100.0, 400.0))
    observations = laboratory ? (SVector(1.5, 0.005, 0.7), SVector(1.7, 0.005, 0.1)) :
        three_dimensional ? (SVector(4500.0, 2500.0, 545.0), SVector(5100.0, 2500.0, -55.0)) :
        (SVector(4500.0, 0.5, 700.0), SVector(5100.0, 0.5, 100.0))
    return (
        physical_size = laboratory ? (2.8, 0.01, 1.2) : (8400.0, three_dimensional ? 5000.0 : 1.0, 1200.0),
        geometry_scale = laboratory ? (1.0, 1.0) : spe11_geometry_scale,
        permeability = permeability, porosity = porosity,
        vertical_permeability_ratio = laboratory ? 1.0 : 0.1,
        exponent = laboratory ? 2.0 : 1.5,
        max_capillary_pressure = laboratory ? 9.5e4 : 3e7,
        entry_pressure = entry_pressure,
        rock_conductivity = (1.9, 1.25, 1.25, 1.25, 0.92, 0.26, 2.0),
        diffusion = (1e-9, laboratory ? 1.6e-5 : 2e-8),
        dispersivity = laboratory ? 0.01 : 10.0,
        rock_density = 2500.0, rock_heat_capacity = 850.0,
        buffer_facies = 2:5,
        well_radius = laboratory ? 9e-4 : 0.15,
        well_positions = laboratory ? (SVector(0.9, 0.005, 0.9), SVector(1.7, 0.005, 0.5)) :
            (SVector(2700.0, 0.5, 900.0), SVector(5100.0, 0.5, 500.0)),
        observation_points = observations, reporting_boxes = boxes,
        datum_depth = laboratory ? 0.0 : 900.0,
        datum_pressure = laboratory ? 1.1e5 : 3e7,
        initialization_time = laboratory ? 0.0 : 1000*spe11_year,
        injection_times = laboratory ? (9000.0, 9000.0) : (25*spe11_year, 25*spe11_year),
        simulation_time = laboratory ? 432000.0 : 1000*spe11_year,
        injection_rate = laboratory ? 1.7e-7 : three_dimensional ? 50.0 : 0.035,
        injection_temperature = laboratory ? 293.15 : 283.15,
        sparse_interval = laboratory ? 600.0 : spe11_year/10,
        report_dims = laboratory ? (280, 1, 120) : three_dimensional ? (168, 100, 120) : (840, 1, 120)
    )
end

spe11_geothermal_temperature(depth) = 343.15 - 0.025*(1200.0 - depth)
spe11_temperature(depth, case) = case == :a ? 293.15 : spe11_geothermal_temperature(depth)

function spe11_reporting_times(case; kind = :dense)
    spec = spe11_case_specification(case)
    if kind == :sparse
        return collect(0.0:spec.sparse_interval:spec.simulation_time)
    elseif kind != :dense
        throw(ArgumentError("Reporting kind must be :dense or :sparse"))
    elseif case == :a
        return collect(0.0:si_unit(:hour):spec.simulation_time)
    elseif case == :b
        return collect(0.0:5*spe11_year:spec.simulation_time)
    else
        return [collect(0:5:50); 75; 100; collect(150:50:500); collect(600:100:1000)].*spe11_year
    end
end
