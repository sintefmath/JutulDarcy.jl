"""
    [rho_brine, c_brine] = pvt_brine_RoweChou1970(T, P, S)

Calculate brine density and/or compressibility using Rowe and Chou"s
(1970) correlation.

# Parameter range for correlation:

P <= 35 MPa
293.15 <= T <= 423.15 K

# Arguments:

  - T: Scalar with temperature value in Kelvin
  - P: Scalar with pressure value in bar
  - S: Salt mass fraction

# Outputs

  - rho_brine:    Scalar with density value in kg/m3
  - c_brine:      Scalar with compressibility value in 1/kPa

"""
function pvt_brine_RoweChou1970(T, P, S; check::Bool = true)
    if check
        # Check if T, P conditions are within range
        if P > 350
            jutul_message("pvt_brine_RoweChou1970", "Pressure $P out of range", color = :yellow)
        end
        if T < 293.15 || T > 423.15
            jutul_message("pvt_brine_RoweChou1970", "Temperature $T out of range 293.15 to 423.15 in pvt_brine_RoweChou1970", color = :yellow)
        end
    end
    bar_to_kgfcm2 = 1.01972
    q = P * bar_to_kgfcm2                        # [kgf/cm^2]
    p_kpa = P * 100                                  # [kPa]

    # Correlation factors
    a1 = 5.916365 - 0.01035794 * T + 0.9270048 * 10^(-5) * T^2 - 1127.522 / T + 100674.1 / T^2
    a2 = 0.520491 * 10^(-2) - 0.10482101 * 10^(-4) * T + 0.8328532 * 10^(-8) * T^2 - 1.1702939 / T + 102.2783 / T^2
    a3 = 0.118547 * 10^(-7) - 0.6599143 * 10^(-10) * T
    a4 = -2.5166 + 0.0111766 * T - 0.170522 * 10^(-4) * T^2
    a5 = 2.84851 - 0.0154305 * T + 0.223982 * 10^(-4) * T^2
    a6 = -0.0014814 + 0.82969 * 10^(-5) * T - 0.12469 * 10^(-7) * T^2
    a7 = 0.0027141 - 0.15391 * 10^(-4) * T + 0.22655 * 10^(-7) * T^2
    a8 = 0.62158 * 10^(-6) - 0.40075 * 10^(-8) * T + 0.65972 * 10^(-11) * T^2

    # Compute density and compressibility
    rho_brine = 1 / ((a1 - q * a2 - q^2 * a3 + a4 * S + a5 * S^2 - q * a6 * S - q * a7 * S^2 - 0.5 * q^2 * a8 * S) * 0.001)                  # [kg/m^3]

    p_ref = 101.325                                   # [kPa]
    r = (p_ref * 0.01) * bar_to_kgfcm2                # [kgf/cm^2]
    rho_ref = 1 / ((a1 - r * a2 - r^2 * a3 + a4 * S + a5 * S^2 - r * a6 * S - r * a7 * S^2 - 0.5 * r^2 * a8 * S) * 0.001)                  # [kgf/m^3]
    c_brine = (rho_brine - rho_ref) / (rho_brine * (p_kpa - p_ref)) # [1/kPa]


    return (rho_brine, c_brine)
end

"""
    mu = viscosity_co2_Fenghour1998(T, rho)

Calculate CO2 viscosity from Vesovic et al., J Phys Chem Ref Data (1990) 
and Fenghour et al., J Phys Chem Ref Data (1998), as described in
Hassanzadeh et al., IJGGC (2008). 


# Arguments:

  - T: Scalar of temperature value in Kelvin
  - rho: Scalar of density value in kg/m^3

# Outputs

  - mu: Dynamic viscosity in Pa*s
"""
function viscosity_co2_Fenghour1998(T, rho; check = true)

    # Check if T, P conditions are within range

    # From Hassanzadeh et al., IJGGC (2008)
    # Constants
    e = [0.235156, -0.491266, 5.211155 * 10^-2,
        5.347906 * 1e-2, -1.537102 * 1e-2]
    f = [5.5934 * 1e-3, 6.1757 * 1e-5, 0.0, 2.6430 * 1e-11]
    g = [0.4071119 * 1e-2, 0.7198037 * 1e-4, 0.2411697 * 1e-16,
        0.2971072 * 1e-22, -0.1627888 * 1e-22]
    exps = [0 1 2 3 4]

    # Functions of scaled T
    T_x = T * (1 / 251.196)
    psi_m = exp(sum(e .* log(T_x) .^ exps))

    # Compute viscosity
    a1 = 1.00697 * T^0.5 / psi_m
    a2 = g[1] * rho
    a3 = g[2] * rho^2
    a4 = g[3] * rho^6 / T_x^3
    a5 = g[4] * rho^8
    a6 = g[5] * rho^8 / T_x
    if T < 313.15
        a7 = sum(f .* rho)
    else
        a7 = 0
    end
    mu = (a1 + a2 + a3 + a4 + a5 + a6 + a7) / 10^6
    return mu
end

"""
    rho_b = pvt_brine_BatzleWang1992(T, P, w_nacl)

Calculate the brine (H2O + NaCl) density based on Batzle & Wang (1992).
These authors used the data of Rowe and Chou (1970), Zarembo & Fedorov
(1975) and Potter & Brown (1977) to expand the P, T validity range. 

# Parameter range for correlation:

P valid from 5 to 100 MPa, T from 20 to 350 C (Adams & Bachu, 2002)

# Arguments

  - T: Temperature value in Kelvin
  - P: Pressure value in bar
  - w_nacl: Salt (NaCl) mass fraction

# Outputs 
  - rho_b: Scalar with brine density in kg/m3
"""
function pvt_brine_BatzleWang1992(T, P, w_nacl; check::Bool = true)

    # Check range
    if check
        if P < 50 || P > 1000
            jutul_message("pvt_brine_BatzleWang1992", "Pressure $P out of measured range 50 to 1000", color = :yellow)
        end
        if T < 273.15 + 20 || T > 273.15 + 350
            jutul_message("pvt_brine_BatzleWang1992", "Temperature $(T-273.15) out of tested range 20 to 350 C", color = :yellow)
        end
        if w_nacl > 0.3 # approx, limit is 320,000 mg/l
            jutul_message("pvt_brine_BatzleWang1992", "Salinity $w_nacl outside of tested range 0 to 0.3", color = :yellow)
        end
    end

    P = 0.1 * P
    T = T - 273.15

    # Water density
    rho_w = 1 + 10^(-6) * (-80 * T - 3.3 * T^2 + 0.00175 * T^3 + 489 * P
                           -
                           2 * T * P + 0.016 * T^2 * P - 1.3 * 10^(-5) * T^3 * P
                           -
                           0.333 * P^2 - 0.002 * T * P^2)                            # [g/cm^3]

    # Brine density
    fact1 = 300 * P - 2400 * P * w_nacl + T * (80 + 3 * T - 3300 * w_nacl - 13 * P + 47 * P * w_nacl)
    rho_b = (rho_w + w_nacl * (0.668 + 0.44 * w_nacl + 10^(-6) * fact1)) * 1000      # [kg/m^3]

    return rho_b
end

"""
    [V_m, rhox, rho] = pvt_co2_RedlichKwong1949(T, P)
    [V_m, rhox, rho] = pvt_co2_RedlichKwong1949(T, P, a_m, b_m)

Calculate CO2 molar volume and density using Redlich and Kwong (1949)
EoS (= RK EoS).

# Parameter range for correlation:

Tested by Spycher et al. (2003) with constant intermolecular attraction
parameters to yield accurate results in the T range ~10 to ~100C and P range up
to 600 bar, for (1) CO2 compressibility factor, (2) CO2 fugacity coefficient and
(3) mutual solubilities of H2O and CO2 in the gas and aqueous phase
(respectively).

# Arguments

  - T: Scalar with temperature value in Kelvin
  - P: Scalar with pressure value in bar

# Optional arguments:

  - a_m: Intermolecular attraction constant (of the mixture) in bar*cm^6*K^0.5/mol^2
  - b_m: Intermolecular repulsion constant (of the mixture) in cm^3/mol

# Outputs

  - V_m: Scalar with molar volume in [cm^3/mol]
  - rhox: Scalar with density in [mol/m^3]
  - rho: Scalar with density in [kg/m^3]
"""
function pvt_co2_RedlichKwong1949(T, P, a_m=7.54 * 10^7 - 4.13 * 10^4 * T, b_m=27.8; check = true)
    # Check if T, P conditions are within range
    if check
        if P > 600
            jutul_message("pvt_co2_RedlichKwong1949", "Pressure $P out of tested range", color = :yellow)
        end
        if P < 0.0
            error("Negative pressure $P")
        end
        if T < 0.0
            error("Negative temperature $T")
        end
        if T < 283.15 || T > 373.15
            jutul_message("pvt_co2_RedlichKwong1949", "Temperature $T out of tested range", color = :yellow)
        end
    end

    # Ideal gas constant
    R = 83.1447                                  # bar*cm^3/mol*K

    # These are for CO2 after Spycher et al. (2003). Infinite dilution of
    # H2O in the gaseous phase is assumed (with y_h2o =0 and y_co2 = 1 in
    # the mixing rules).

    # Molar volume 
    #   Make sure not to leave gaps between signs and each coefficient
    #   when using the "roots" fcn. Also, more than one value below the 
    #   critical point is returned. The gas phase volume is always given by the 
    #   maximum root (see Appx B.2, Spycher et al., 2003).
    c3 = 1
    c2 = (R * T / P)
    c1 = ((R * T * b_m) / P - a_m / (P * T^0.5) + b_m^2)
    c0 = a_m * b_m / (P * T^0.5)
    # r = roots([c3 -c2 -c1 -c0]);
    r = roots(Polynomial([-c0, -c1, -c2, c3])) # TODO: Check this
    # V_m = r(abs(imag(r))<eps);                              # [cm^3/mol]
    V_m = Float64[]
    for root in r
        if abs(imag(root)) < 1e-12
            push!(V_m, root)
        end
    end
    if length(V_m) > 1
        epsval = 1e-12
        V_gas = maximum(V_m)                                 # V gas phase  
        V_liq = minimum(V_m)                                 # V of liquid phase
        w1 = P * (V_gas - V_liq)
        w21 = R * T * log((V_gas - b_m) / (V_liq - b_m))
        w22 = a_m / (b_m * T^0.5) * log(((V_gas + b_m) * V_liq) / ((V_liq + b_m) * V_gas))
        w2 = w21 + w22
        if w2 - w1 > 0 + epsval
            V_m = V_gas
        elseif w2 - w1 < 0 - epsval
            V_m = V_liq
        else
            V_m = V_gas
        end
    else
        V_m = only(V_m)
    end

    # Molar density
    rhox = (1 / V_m) * 10^6                                    # [mol/m^3]

    # Mass density
    Mw_co2 = 44.0095 / 10^3                                  # [kg/mol]
    rho = rhox * Mw_co2                                   # [kg/m^3]
    return (V_m, rhox, rho)
end

"""
    gamma_co2 = activity_co2_DS2003(T, P, m_io)

Calculate a CO2 pseudo activity coefficient based on a virial expansion
of excess Gibbs energy.

# Arguments
  - T: Scalar with temperature value in Kelvin
  - P: Scalar with pressure value in bar
  - m_io: Vector where each entry corresponds to the
    molality of a particular ion in the initial brine solution. The
    order is as follows:
    [ Na(+),   K(+),  Ca(2+), Mg(2+), Cl(-), SO4(2-)]

# Outputs 
  - V_m: Scalar with molar volume in [cm^3/mol]
  - rhox: Scalar with density in [mol/m^3]
  - rho: Scalar with density in [kg/m^3]
"""
function activity_co2_DS2003(T, P, m_io; check = true)
    if check
        # Check if T, P conditions are within range
        if P > 2000
            jutul_message("activity_co2_DS2003", "Pressure out of tested range", color = :yellow)
        end
        if T < 273 || T > 533
            jutul_message("activity_co2_DS2003", "Temperature out of tested range", color = :yellow)
        end
        if maximum(m_io) > 4.3
            jutul_message("activity_co2_DS2003", "Ionic strength out of tested range", color = :yellow)
        end
    end

    # Interaction parameter constants
    c1_co2_C = -0.411370585
    c2_co2_C = 6.07632013 * 10^(-4)
    c3_co2_C = 97.5347708
    c4_co2_C = -0.0237622469
    c5_co2_C = 0.0170656236
    c6_co2_C = 1.41335834 * 10^(-5)

    c1_co2_C_A = 3.36389723 * 10^(-4)
    c2_co2_C_A = -1.98298980 * 10^(-5)
    c3_co2_C_A = 0
    c4_co2_C_A = 2.122220830 * 10^(-3)
    c5_co2_C_A = -5.24873303 * 10^(-3)
    c6_co2_C_A = 0

    # Compute interaction parameters
    lam_co2_C = c1_co2_C + c2_co2_C * T + c3_co2_C / T +
                c4_co2_C * P / T + c5_co2_C * P / (630 - T) +
                c6_co2_C * T * log(P)

    zet_co2_C_A = c1_co2_C_A + c2_co2_C_A * T + c3_co2_C_A / T +
                  c4_co2_C_A * P / T + c5_co2_C_A * P / (630 - T) +
                  c6_co2_C_A * T * log(P)

    # Compute activity coefficient
    gamma_co2 = exp(2 * lam_co2_C * (sum(m_io[1:2]) + 2 * sum(m_io[3:4])) + zet_co2_C_A * m_io[5] * (sum(m_io[1:4])) - 0.07 * m_io[6])
    return gamma_co2
end

function convert_salinity_to_mole_fractions(salinity)
    mw_h2o = 0.01801528
    moles_nacal = 1e3*salinity*1.0
    moles_water = 1000.0/mw_h2o
    mf_nacl = moles_nacal/(moles_nacal + moles_water)
    return mf_nacl
end

function compute_salinity(inp_mole_fractions=Float64[], names=String[])
    @assert sum(inp_mole_fractions) < 0.9 "Should be mass fractions without water"
    saltSpecies = ("NaCl", "KCl", "CaSO4", "CaCl2", "MgSO4", "MgCl2")
    nu = (2, 2, 2, 3, 2, 3)
    ions = ("Na_+", "K_+", "Ca_2+", "Mg_2+", "Cl_-", "SO4_2-")
    #         [ NaCl,     KCl,    CaSO4,  CaCl2,  MgSO4,  MgCl2]
    Mw_salt = (58.44277, 74.5513, 136.1406, 110.984, 120.3676, 95.211) ./ 1e3       # [kg/mol)
    Mw_h2o = 18.01528 / 1e3
    #         [ Na(+),   K(+),  Ca(2+), Mg(2+), Cl(-), SO4(2-)]
    Mw_io = (22.98977, 39.0983, 40.078, 24.305, 35.453, 96.0626) ./ 1e3           # [")


    mole_fractions = zeros(6)
    for i in eachindex(inp_mole_fractions, names)
        name = names[i]
        ix = findfirst(isequal(name), saltSpecies)
        @assert !isnothing(ix) "$name not in $saltSpecies"
        mole_fractions[ix] = inp_mole_fractions[i]
    end

    mass_salts = mole_fractions .* Mw_salt
    mass_salt = sum(mass_salts)
    h2o_mass = (1.0 - sum(mole_fractions)) * Mw_h2o

    w_salt = mass_salt ./ (h2o_mass + mass_salt)
    # m_salt is molality of species k in solution (mol k / kg solvent)
    m_salt = mole_fractions ./ h2o_mass
    m_io = (m_salt[1], m_salt[2], m_salt[3] + m_salt[4], m_salt[5] + m_salt[6], sum(m_salt[1:2]) + 2 * (m_salt[4] + m_salt[6]), m_salt[3] + m_salt[5])
    return (m_salt, m_io, nu, Mw_salt, Mw_io, w_salt)
end


"""
    mu_b_co2 = viscosity_brine_co2_mixture_IC2012(T, P, m_nacl, w_co2)

Calculate the dynamic viscosity of a solution of H2O + NaCl (brine) with
dissolved CO2.

# Parameter range for correlation:

For pure water + CO2, the model is based on experimental data by Kumagai
et al. (1998), which is for p up to 400 bar and T up to 50 C, and Bando 
et al. (2004), which is valid in 30-60C and 10-20MPa.
The model of Mao & Duan (2009) for brine viscosity reaches 623K, 1000 bar
and high ionic strength. However, the model used to determine the viscosity 
when co2 dissolves in the brine (Islam & Carlson, 2012) is based on
experimental data by Bando et al. (2004) and Fleury and Deschamps (2008),
who provided experimental data up to P = 200 bar, T extrapolated to
100 C, and maximum salinity of 2.7M.

# Arguments

  - T: Scalar with temperature value in Kelvin
  - P: Scalar with pressure value in bar
  - m_nacl: Salt molality (NaCl) in mol/kg solvent
  - w_co2: Mass fraction of CO2 in the aqueous solution (i.e. brine)

# Outputs

  - mu_b_co2: Scalar with dynamic viscosity in Pa*s

"""
function viscosity_brine_co2_mixture_IC2012(T, P, m_nacl, w_co2; check=true)
    # Check range
    if check
        if m_nacl == 0
            if P > 400
                jutul_message("viscosity_brine_co2_mixture_IC2012", "Pressure out of measured range", color = :yellow)
            end
            if T < 273.15 || T > 333.15
                jutul_message("viscosity_brine_co2_mixture_IC2012", "Temperature out of tested range", color = :yellow)
            end
        else
            if P > 200
                jutul_message("viscosity_brine_co2_mixture_IC2012", "Pressure out of measured range", color = :yellow)
            end
            if T < 273.15 + 35 || T > 373.15
                jutul_message("viscosity_brine_co2_mixture_IC2012", "Temperature out of tested range", color = :yellow)
            end
            if m_nacl > 3.1 # m approx, limit of experimental data is 2.738M
                jutul_message("viscosity_brine_co2_mixture_IC2012", "Salinity of tested range", color = :yellow)
            end
        end
    end
    # Units
    P_mpa = P / 10.0

    # Pure water density (see Islam & Carlson, 2012)
    a = 1.34136579 * 10^2
    b = [-4.07743800 * 10^3, 1.63192756 * 10^4, 1.37091355 * 10^3]
    c = [-5.56126409 * 10^-3, -1.07149234 * 10^-2, -5.46294495 * 10^-4]
    d = [4.45861703 * 10^-1, -4.51029739 * 10^-4]
    cp = [1, 2]

    rho_h2o = a + sum(b .* 10 .^ (c .* T)) + sum(d .* P_mpa .^ cp)                     # [kg/m3]
    rho_h2o = rho_h2o / 10^3                                                   # [g/cm^3]

    # Pure water viscosity (Mao & Duan, 2009)
    d = [
        0.28853170 * 10^7, -0.11072577 * 10^5, -0.90834095 * 10,
        0.30925651 * 10^-1, -0.27407100 * 10^-4, -0.19283851 * 10^7,
        0.56216046 * 10^4, 0.13827250 * 10^2, -0.47609523 * 10^-1,
        0.35545041 * 10^-4
    ]
    coefs1 = (1:5) .- 3
    coefs2 = (6:10) .- 8

    mu_h2o = exp(sum(d[1:5] .* T .^ coefs1) + sum(d[6:10] .* rho_h2o .* T .^ coefs2))

    if m_nacl == 0
        # Viscosity of H2O + CO2 (Islam & Carlson, 2012)
        # a = [7.632609119*1e2 -9.46077673946*1e3]; # Reported by IC2012
        a = (1.632609119 * 1e3, -9.46077673946 * 1e2)   # Actual fit to plots in paper + data by Bando et al. (2004)
        b = (-1.047187396332 * 1e4, 3.68325597 * 1e1)
        c1 = (1, 2)
        c2 = (0, 1)
        mu_r = 1 + sum(a .* w_co2 .^ c1) / sum(b .* T .^ c2)
        mu_b_co2 = mu_r * mu_h2o
    else
        # Brine viscosity (H2O + NaCl) ( Mao & Duan, 2009)
        A = -0.21319213 + 0.13651589 * 10^(-2) * T - 0.12191756 * 10^(-5) * T^2
        B = 0.69161945 * 10^-1 - 0.27292263 * 10^(-3) * T + 0.20852448 * 10^(-6) * T^2
        C = -0.25988855 * 10^-2 + 0.77989227 * 10^(-5) * T
        mu_b = exp(A * m_nacl + B * m_nacl^2 + C * m_nacl^3) * mu_h2o
        # Brine viscosity with dissolved CO2 (H2O + NaCl + CO2) (Islam & Carlson,
        # 2012)
        mu_b_co2 = mu_b * (1 + 4.65 * w_co2^1.0134)
    end
    # Data from Bando et al., J Chem Eng Data (2004) for test
    # xcs = P_mpa/(36.1*P_mpa +3.87*T-1097.1+w_salt*(196*P_mpa + 26.9*T-8810))
    # dat = [303.15 100 0.89e-3; 303.15 200 0.88e-3; 313.15 100 0.7e-3; 
    #        313.15 200 0.71e-3; 323.15 100 0.57e-3; 323.15 200 0.57e-3;
    #        333.15 100 0.47e-3; 333.15 200 0.48e-3];
    return mu_b_co2
end

"""
    viscMixture = viscosity_gas_mixture_Davidson1993(x, M, mu)

Calculate the viscosity of a gas mixture following Davidson (1993). In
principle, valid for any range within which the individual components'
viscosities are valid.

Arguments:

  - x: Mole fraction of each component
  - M: Molar mass of each component
  - mu: Viscosity of each component in user chosen units

Each input should be a Float64 Vector of length n where n is the total
number of components

# Outputs

  - viscMixture: Scalar viscosity of the mixture in same units as mu
"""
function viscosity_gas_mixture_Davidson1993(x, M, mu)

    # Momentum fractions
    num = x .* sqrt(M)
    den = sum(num)
    y = num / den

    # Mass fraction
    M_mixt = sum(x .* M)
    w = x .* (M ./ M_mixt)

    # Efficiency of momentum transfer
    E = 2 .* sqrt(w)' .* sqrt(w) ./ (w' + w)
    for (i, eval) in enumerate(E)
        if isnan(eval)
            E[i] = 0
        end
    end
    # E(isnan(E)) = 0;                   # fix for when some component is absent.

    # Fluidity
    A = 0.375  # Best fit (N = 164 mixtures) with RMS = 1.28# (Davidson, 1993)
    f = sum(sum(y' .* y ./ (sqrt(mu)' .* sqrt(mu)) .* E .^ A))

    # Viscosity
    viscMixture = 1 / f

    # (Slow way of computing E and F)
    # N = numel(x);
    # f = 0;
    # for i = 1:N
    #     for j = 1:N
    #         E = 2*sqrt(w(i))*sqrt(w(j))/(w(i) + w(j));
    #         f = f + (y(i)*y(j)/(sqrt(mu(i))*sqrt(mu(j)))*E^A);
    #     end
    # end
    return viscMixture
end

"""
    props = compute_co2_brine_props(p_pascal, T_K, salt_mole_fractions = Float64[], salt_names = String[];
        check=true,
        iterate=false,
        maxits=15,
        ionized=false
    )
    props = compute_co2_brine_props(200e5, 273.15 + 30.0)

Get pure phase properties (density, viscosity) and equilibrium constants for
brine and water. This functions and the functions used in this function are
heavily based on comparable code in the the co2lab-mit MRST module developed by
Lluis Salo.

The salt mole fractions are of the total brine (i.e. including H2O in the
calculation) and can any subset of the following names provided:

    ("NaCl", "KCl", "CaSO4", "CaCl2", "MgSO4", "MgCl2")

"""
function compute_co2_brine_props(p_pascal, T_K, salt_mole_fractions = Float64[], salt_names = String[];
        check = true,
        iterate = false,
        verbose = false,
        maxits = 15,
        ionized = false
    )
    # TODO: Salt var name
    saltVar = ionized
    m_salt, m_io, nu, Mw_salt, Mw_io, w_salt = compute_salinity(salt_mole_fractions, salt_names)
    p = p_pascal / 1e5
    P = p
    T_C = T_K - 273.15
    T_K_s = 273.15 + 15.56
    P_bar_s = 1.013529

    if sum(m_salt[2:end]) ≈ 0.0
        pvt_brine = pvt_brine_RoweChou1970
    else
        pvt_brine = pvt_brine_BatzleWang1992
    end
    rho_b_fcn(arg...) = pvt_brine(arg...; check = check)

    rho_co2_fcn(arg...) = pvt_co2_RedlichKwong1949(arg...; check = check)
    rho_gas_fcn(Mw_gas, V) = Mw_gas / V
    mu_aq_fcn(arg...) = viscosity_brine_co2_mixture_IC2012(arg...; check = check)
    mu_gas_fcn(arg...) = viscosity_gas_mixture_Davidson1993(arg...; check = check)
    mu_co2_fcn(arg...) = viscosity_co2_Fenghour1998(arg...; check = check)


    if check
        if p > 600 && p < 700
            jutul_message("compute_co2_brine_props", "Input pressure contains values above 600 bar. Model may not be accurate.", color = :yellow)
        elseif p > 700
            jutul_message("compute_co2_brine_props", "Input pressure contains values significantly above 600 bar. Formulation possibly not accurate.", color = :red)
        end
        if T_C > 100 && T_C <= 125
            jutul_message("compute_co2_brine_props", "Input temperature above 100 C. Model may not be accurate.", color = :yellow)
        elseif T_C > 125
            jutul_message("compute_co2_brine_props", "Input temperature significantly above 100 C. This formulation is not accurate.", color = :red)
        end
        if maximum(m_salt) > 4 && maximum(m_salt) < 6
            jutul_message("compute_co2_brine_props", "Input salinity above 4 m. Model may not be accurate.", color = :yellow)
        elseif maximum(m_salt) > 6
            jutul_message("compute_co2_brine_props", "Input salinity above halite saturation. This formulation is not accurate.", color = :red)
        end
    end

    R = 83.1447                                          # bar*cm^3/mol*K
    Mw_h2o = 18.01528 / 10^3                                    # [kg/mol]
    Mw_co2 = 44.0095 / 10^3                                     #   "

    # Compute mole fractions
    m_h2o = 1 / Mw_h2o                                         # mol H20/ kg H2O
    vm = nu .* m_salt
    if saltVar
        x_salt = (vm) ./ (sum(vm) + m_h2o)    # [-] mole fractions (total), fully ionized salt
        x_salts = (vm) ./ sum(vm)              # [-] (of salts)
    else
        x_salt = m_salt ./ (sum(m_salt) + m_h2o)              # [-] Not fully ionized (Hassanzadeh et al., 2008)
        x_salts = m_salt ./ sum(m_salt)                          # [-] 
    end

    x_h2o = 1 - sum(x_salt)                                  # [-]
    Mw_salts = sum(Mw_salt .* x_salts)                           # [kg/mol] avg

    # Compute mixture Mw (Mw_brine) in kg/mol
    Mw_b = sum(Mw_salt .* x_salt) + Mw_h2o * x_h2o                 # [kg/mol]

    # Water and CO2 Redlich-Kwong molecular interaction params (Spycher et al., 2003) 
    # The parameters involving H20 were calculated assuming infinite dilution
    # of water in the gas phase, so they should be recalculated if this
    # assumption is relaxed with this solubility model. In this case, a_h20
    # would also be needed (not reported by Spycher et al., 2003, due to this
    # assumption).
    a_co2 = 7.54 * 10^7 - 4.13 * 10^4 * T_K                    # bar*cm^6*K^0.5/mol^2
    a_h2o_co2 = 7.89 * 10^7                                    # "
    b_co2 = 27.8                                         # cm^3/mol
    b_h2o = 18.18                                        # "
    a_h2o = 0                                            # not used here

    # Average partial molar volumes of each pure condensed component, they are
    # assumed constant in the p interval of interest (Spycher et al., 2003)
    Vp_co2 = 32.6                                         # cm^3/mol (both (g) and (l))
    Vp_h2o = 18.1                                         # "

    # True equilibrium ctnts. at p0=1bar (/kappa params; Spycher et al., 2003)
    p0 = 1                                            # bar
    kap0_co2_g = 10^(1.189 + 1.304 * 10^(-2) * T_C - 5.446 * 10^(-5) * T_C^2)
    kap0_co2_l = 10^(1.169 + 1.368 * 10^(-2) * T_C - 5.380 * 10^(-5) * T_C^2)
    kap0_h2o = 10^(-2.209 + 3.097 * 10^(-2) * T_C - 1.098 * 10^(-4) * T_C^2 + 2.048 * 10^(-7) * T_C^3)

    # Partial molar volume of CO2 (Garcia, 2001)
    V_phi = 37.51 - 9.585 * 10^(-2) * T_C + 8.740 * 10^(-4) * T_C^2 - 5.044 * 10^(-7) * T_C^3                          # [cm^3/mol]
    V_phi = V_phi / 10^6                                     # [m^3/mol]

    # CO2 and brine molar densities at standard conditions
    Y0_h2o = 0
    Y_h2o_out = 0
    Y0_co2 = 1
    Y_co2_out = 1

    # Compute
    idt2 = 1
    idp = 0
    it = 0

    Y_co2 = X_co2 = rho_brine = rho_gas = rho_co2 = NaN
    while it <= maxits
        # Mixing rules (Prausnitz et al., 1986)
        a_m = Y0_h2o^2 * a_h2o + 2 * Y0_co2 * Y0_h2o * a_h2o_co2 + Y0_co2^2 * a_co2
        b_m = Y0_co2 * b_co2 + Y0_h2o * b_h2o

        # 1. Gas molar volume (Redlich and Kwong (1949) EoS, = RK EoS)
        V, _, rho_co2 = pvt_co2_RedlichKwong1949(T_K, P, a_m, b_m, check = check)              # [cm^3/mol, ~, kg/m^3]

        # 2. Gas compressibility factor
        Z = P * V / (R * T_K)                                              # [-]

        # 3. Fugacity coefficients (Spycher et al., 2003)
        #   3.1 CO2
        fa = log(V / (V - b_m))
        fb = b_co2 / (V - b_m)
        fc = (2 * (Y0_co2 * a_co2 + Y0_h2o * a_h2o_co2) / (R * T_K^1.5 * b_m)) * log((V + b_m) / V)
        fd = a_m * b_co2 / (R * T_K^1.5 * b_m^2) * (log((V + b_m) / V) - b_m / (V + b_m))
        fe = log(Z)
        phi_co2 = exp(fa + fb - fc + fd - fe)                        # [-]

        #   3.2 Water
        fbw = b_h2o / (V - b_m)
        fcw = (2 * (Y0_co2 * a_h2o_co2 + Y0_h2o * a_h2o) / (R * T_K^1.5 * b_m)) * log((V + b_m) / V)
        fdw = a_m * b_h2o / (R * T_K^1.5 * b_m^2) * (log((V + b_m) / V) - b_m / (V + b_m)) # typo in Hassanzadeh et al., 2008
        phi_h2o = exp(fa + fbw - fcw + fdw - fe)                     # [-]

        # 4. CO2 molality in pure H2O at P, T conditions (m0_co2)
        if T_C < 31 && V < 94                                         # Spycher et al., 2003
            kap0_co2 = kap0_co2_l
        else
            kap0_co2 = kap0_co2_g
        end
        B = phi_co2 * P / (m_h2o * kap0_co2) * exp(-(P - p0) * Vp_co2 / (R * T_K))
        A = kap0_h2o / (phi_h2o * P) * exp((P - p0) * Vp_h2o / (R * T_K))
        # y_i = mole fraction of component i in the gaseous phase
        y_h2o = (1 - B) / (1 / A - B)                                        # [-] Spycher et al., 2003
        # x_i = mole fraction of component i in the aqueous phase
        x_co2 = B * (1 - y_h2o)                                          # [-]
        # Check values within range
        if x_co2 < 0 && check
            @warn "x_co2 = $x_co2. Set to 0"
            x_co2 = 0
        elseif x_co2 > 1 && check
            @warn "x_co2 = $x_co2. Set to 1"
            x_co2 = 1
        end
        if y_h2o < 0 && check
            @warn "y_h2o = $y_h2o. Set to 0"
            y_h2o = 0
        elseif y_h2o > 1 && check
            @warn "y_h2o = $y_h2o. Set to 1"
            y_h2o = 1
        end
        x_h2o = 1 - x_co2                                            # [-]
        # molality
        m0_co2 = m_h2o * x_co2 / x_h2o                                   # [m]

        # 5. CO2 molality in saline solution at P, T conditions (m_co2)
        #   5.1 CO2 activity coefficient (Duan and Sun, 2003, as in Spycher
        #   & Pruess, 2005; typos in Hassanzadeh et al., 2008)
        gamma_co2 = activity_co2_DS2003(T_K, P, m_io, check = check)

        #   5.2 CO2 molality
        m_co2 = m0_co2 / gamma_co2

        # 6. H2O and CO2 mole fractions and equilibrium ratios in aqueous
        # saline solution with dissolved CO2.
        X_co2 = m_co2 / (m_co2 + m_h2o + sum(vm))              # [-] mole fraction
        if saltVar
            X_salt = sum(vm) / (m_co2 + m_h2o + sum(vm))       # [-] fully ionized
        else
            X_salt = sum(x_salt)                                     # [-] Not fully ionized
        end
        X_solv = 1 - X_co2                                        # [-]
        X_h2o = X_solv - X_salt                                  # [-] fully ionized
        Y_h2o = A * X_h2o                                          # [-] fully ionized
        Y_co2 = 1 - Y_h2o                                          # [-]
        K_co2 = Y_co2 / X_co2                                      # [-] equilibrium rat.
        K_h2o = Y_h2o / X_h2o                                      # [-]
        V, _, rho_co2 = rho_co2_fcn(T_K, P, a_m, b_m)
        rho_brine, = rho_b_fcn(T_K, P, w_salt)
        rho_gas = rho_gas_fcn(Mw_co2, V / 10^6)
        mu_brine = mu_aq_fcn(T_K, P, 0, 0)
        mu_gas = mu_co2_fcn(T_K, rho_co2)

        done = false
        if iterate
            res = abs(Y_h2o - Y0_h2o)
            Y0_h2o = Y_h2o
            Y0_co2 = Y_co2
            if res < 1e-6
                verbose && println("Iterations at p=$P bar: $it")
                done = true
            elseif it == maxits
                error("Iterative loop did not converge")
            end
            it = it + 1
        else
            done = true
        end
        if done
            out = Dict{Symbol,Any}()
            out[:density] = @SVector [rho_brine, rho_gas]
            out[:viscosity] = @SVector [mu_brine, mu_gas]
            out[:K] = @SVector [K_h2o, K_co2]
            return out
        end
    end
    error("Something went wrong, this should not be reached.")
end
