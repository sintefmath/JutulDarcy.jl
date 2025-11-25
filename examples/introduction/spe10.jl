
using HYPRE, JutulDarcy, GLMakie
case = JutulDarcy.SPE10.setup_case(layers = 1:85)
ws, states = simulate_reservoir(case, info_level = 1)

