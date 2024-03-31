# Switch to local
import Pkg

cd()
Pkg.activate("Dynamical_Core_Julia")
cd("Dynamical_Core_Julia")

# import local
using dynamical_core

include("__main.jl")

# execute
Shallow_Water_Main()