__precompile__(true)
module dynamical_core

using Revise
using LinearAlgebra
using FFTW
using Statistics 
using JLD2
using MAT
using Random
import PyPlot

"""
TODO : Seperate SW and generalize method
"""

include("shared/Atmo_Data.jl")
include("shared/Dyn_Data.jl")
include("shared/Gauss_And_Legendre.jl")
include("shared/Spectral_Spherical_Mesh.jl")
include("shared/Vert_Coordinate.jl")
include("shared/Time_Integrator.jl")
include("shared/Semi_Implicit.jl")
include("shared/Press_And_Geopot.jl")

end
