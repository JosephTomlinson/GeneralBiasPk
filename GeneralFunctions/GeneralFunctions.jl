module GeneralFunctions
using Dierckx
using TwoFAST

export BiasDict13, BiasDictlhd, BiasDict22
export extend_func, eval_extend_func
export Plhd, Pllhd, Plhdmultipole
export ξ
export best_q, henry_q, show_best_q

include("BiasDictConstructors.jl")
include("ExtendedFunctions.jl")
include("Plhd.jl")
include("xi.jl")
include("Bestq.jl")
end