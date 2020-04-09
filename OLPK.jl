__precompile__()
module OLPK

using Dierckx
using TwoFAST
using OffsetArrays
using Memoization #Need to go through and memoize many fftlog functions, in progress (Seems to throw an error look into)

const N = 16*1024
const kmino = exp(-25)
const kmaxo = exp(25)
const rmini = exp(-11)#exp(-14.555)
const rmaxi = exp(11)#exp(11)
const r0o = 1/kmaxo
const r0i = 1/rmaxi

export BiasDict13, BiasDictlhd, BiasDict22
export extend_func, eval_extend_func
export Plhd, Pllhd, Plhdmultipole
export ξ
export best_q, henry_q, show_best_q
export Legendre

include("GeneralFunctions/BiasDictConstructors.jl")
include("GeneralFunctions/ExtendedFunctions.jl")
include("GeneralFunctions/Plhd.jl")
include("GeneralFunctions/xi.jl")
include("GeneralFunctions/Bestq.jl")
include("GeneralFunctions/Legendre.jl")

export generateCln
export M

include("P13/Cln.jl")
include("P13/M13.jl")

export generateMnli, generateMnlismall, generateMnliFast

include("P22/Mnli.jl")
include("P22/MnliFast.jl")

export GrowthPowerDict13, P_13ln, I1, I2, I3, I4, I5, P13ggs, P13lggs, P13ggsmultipole

include("P13/P13.jl")

export P22lggs, P22ggsmultipole, intMnlitesting

include("P22/P22.jl")

export BiasDictAll

function BiasDictAll(;b1, bη, bK², bδη, bη², bKKpara, btd, bδΠ2para, bηΠ2para, bΠ2Kpara, bΠ2para, bΠ3para,
    P0_ϵ, b∇2δ, β∇2v, β∂2parav, P2_ϵ, P2_ϵε_η, b2)
    bd13 = BiasDict13(b1=b1, bη=bη, bK²=bK², bδη=bδη, bη²=bη², bKKpara=bKKpara, btd=btd, bδΠ2para=bδΠ2para, 
        bηΠ2para=bηΠ2para, bΠ2Kpara=bΠ2Kpara, bΠ2para=bΠ2para, bΠ3para=bΠ3para)
    bd22 = BiasDict22(b1=b1, bη=bη, b2=b2, bK²=bK², bδη=bδη, bη²=bη², bKKpara=bKKpara, bΠ2para=bΠ2para)
    bdlhd = BiasDictlhd(b1=b1, bη=bη, P0_ϵ=P0_ϵ, b∇2δ=b∇2δ, β∇2v=β∇2v, β∂2parav=β∂2parav, P2_ϵ=P2_ϵ, P2_ϵε_η=P2_ϵε_η)
    
    return merge(bd13, bd22, bdlhd)
end

export PNLOggs

function PNLOggs(Pk, μs; GeneralBiasDict, f)
    k, P22 = P22ggsmultipole(Pk, μs, BiasOperatorDict22=GeneralBiasDict, f=f)
    _, P13 = P13ggs(Pk, μs, BiasOperatorDict13=GeneralBiasDict, f=f)
    Pl = Plhd(Pk, k, μs, BiasOperatorDictlhd=GeneralBiasDict, f=f)
    PNLO = @. Pl + P22 + 2P13
    
    return k, PNLO
end

export PNLOlggs

function PNLOlggs(Pk, l; GeneralBiasDict, f)
    r, Mnli = generateMnliFast(Pk, GeneralBiasDict, f)
    k, P22l = P22lggs(Pk, l, BiasOperatorDict22=GeneralBiasDict, f=f, Mnli=Mnli, r=r)
    
    Cln = generateCln(GeneralBiasDict, f)
    _, P13l = P13lggs(Pk, l, BiasOperatorDict13=GeneralBiasDict, f=f, Cln=Cln)
    
    Pll = Pllhd(Pk, k, l, BiasOperatorDictlhd=GeneralBiasDict, f=f)
    
    PNLOl = @. Pll + P22l + 2P13l
    
    return k, PNLOl
end


end