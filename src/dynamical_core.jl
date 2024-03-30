module dynamical_core
"""
TODO : Seperate SW and generalize method
"""

using FFTW
using JLD2

################################################################################
# Shared
include("shared/Atmo_Data.jl")
include("shared/Dyn_Data.jl")
include("shared/Gauss_And_Legendre.jl")
include("shared/Spectral_Spherical_Mesh.jl")
include("shared/Vert_Coordinate.jl")
include("shared/Time_Integrator.jl")
include("shared/Semi_Implicit.jl")
include("shared/Press_And_Geopot.jl")

################################################################################
# Shallow Water System
include("shallow_water/Output_Manager.jl")
include("shallow_water/Shallow_Water_Physics.jl")
include("shallow_water/Shallow_Water_Dynamics.jl")
include("shallow_water/Initial_Condition_Controller.jl")

################################################################################
# ?

end
