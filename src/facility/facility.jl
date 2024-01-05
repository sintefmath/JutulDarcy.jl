"""
(Absolute) Minimum well rate for a well that is not disabled.
"""
const MIN_ACTIVE_WELL_RATE = 1e-20
"""
(Absolute) Minimum initial rate for wells when controls are updated.
"""
const MIN_INITIAL_WELL_RATE = 1e-12
"""
Well variables - entities that we have exactly one of per well (and usually relates to the surface connection)
"""

include("types.jl")
include("wells/wells.jl")
include("controls.jl")
include("wellgroups.jl")
include("cross_terms.jl")
