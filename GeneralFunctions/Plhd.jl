function Plhd(Pk,ks,μs; BiasOperatorDictlhd, f)
    plhd = zeros(Float64,length(ks),length(μs))
    b1 = BiasOperatorDictlhd["b1"]
    bη = BiasOperatorDictlhd["bη"]
    P0_ϵ = BiasOperatorDictlhd["P0ϵ"]
    b∇2δ = BiasOperatorDictlhd["b∇2δ"]
    β∇2v = BiasOperatorDictlhd["β∇2v"]
    β∂2parav = BiasOperatorDictlhd["β∂2∥v"]
    P2_ϵ = BiasOperatorDictlhd["P2ϵ"]
    P2_ϵε_η = BiasOperatorDictlhd["P2ϵεη"]
    for i = 1:length(μs)
        @. plhd[:,i] = (b1 - bη*f*μs[i]^2)^2*Pk(ks) + P0_ϵ -
        2*(b1*b∇2δ - μs[i]^2*f*bη*(b∇2δ + b1*β∇2v + b1*β∂2parav*μs[i]^2) +
            μs[i]^4*f^2*bη^2*(β∇2v + β∂2parav*μs[i]^2))*ks^2*Pk(ks) +
        ks^2*P2_ϵ + μs[i]^2*ks^2*bη*P2_ϵε_η
    end
    return plhd
end

function Pllhd(Pk, ks, l; BiasOperatorDictlhd, f)
    b1 = BiasOperatorDictlhd["b1"]
    bη = BiasOperatorDictlhd["bη"]
    P0_ϵ = BiasOperatorDictlhd["P0ϵ"]
    b∇2δ = BiasOperatorDictlhd["b∇2δ"]
    β∇2v = BiasOperatorDictlhd["β∇2v"]
    β∂2parav = BiasOperatorDictlhd["β∂2∥v"]
    P2_ϵ = BiasOperatorDictlhd["P2ϵ"]
    P2_ϵε_η = BiasOperatorDictlhd["P2ϵεη"]
    Pkvals = Pk.(ks)
    if l == 0
        @. b1^2*Pkvals - (2/15)b1*Pkvals*(f*bη*(-5β∇2v*ks^2 - 3β∂2parav*ks^2 + 5) + 15b∇2δ*ks^2) -
            (1/35)f^2*Pkvals*bη^2*(14β∇2v*ks^2 + 10β∂2parav*ks^2 - 7) + 
            (1/3)ks^2*bη*(2b∇2δ*f*Pkvals + P2_ϵε_η) + ks^2*P2_ϵ + P0_ϵ
    elseif l == 2
        @. (2/21)bη*(2*f*Pkvals*(f*bη*(3-ks^2*(6β∇2v+5β∂2parav))+b1*(ks^2*(7β∇2v + 6β∂2parav) - 7)) +
            7ks^2*(2b∇2δ*f*Pkvals + P2_ϵε_η))
    elseif l == 4
        @. -(8/385)f*Pkvals*bη*(f*bη*(22β∇2v*ks^2 + 30β∂2parav*ks^2 - 11) - 22b1*β∂2parav*ks^2)
    elseif l == 6
        @. -(32/231)β∂2parav*f^2*ks^2*Pkvals*bη^2
    else
        fill!(similar(ks), 0)
    end
end

function Plhdmultipole(Pk, ks, μs; BiasOperatorDictlhd, f)
    L(n,x) = Legendre(n,x)
    Pl = l -> Pllhd(Pk, ks, l, BiasOperatorDictlhd=BiasOperatorDictlhd, f=f)
    P0 = Pl(0)
    P2 = Pl(2)
    P4 = Pl(4)
    P6 = Pl(6)
    plhdmp = Array{Float64}(undef, length(ks), length(μs))
    for i = 1:length(μs)
        @. plhdmp[:,i] = L(0,μs[i])*P0 + L(2,μs[i])*P2 + L(4,μs[i])*P4 + L(6,μs[i])*P6
    end
    return plhdmp
end