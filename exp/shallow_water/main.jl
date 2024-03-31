# Switch to local
import Pkg

cd()
Pkg.activate("dynamical_core")
cd("dynamical_core")

# import local
using dynamical_core

include("__main.jl")

# execute
Shallow_Water_Main()