function mcP_Lv(Pk) #σᵥ² needed for P13 shortcut, no bias
    f = x -> exp(x)*Pk(exp(x))/(2π^2)
    quadgk(f,log(kmino), log(kmaxo))
end

function P13BiasTest(Pk)
    k, i1 = I1(Pk)
    _, i2 = I2(Pk)
    _, i3 = I3(Pk)
    p13 = @. Pk(k)*3*((2/63)*i1 + (1/42)*i2 - (1/18)*i3 - (1/18)*k^2*mcP_Lv(Pk)[1])
    return k , 2 .*p13
end

function P13BiasTest2(Pk)
    local13BiasDict = BiasDict13(b1=1, bη=0, bK²=0, bδη=0, bη²=0, bKKpara=0, btd=0,
                                 bδΠ2para=0, bηΠ2para=0, bΠ2Kpara=0, bΠ2para=0, bΠ3para=0)
    locallhdBiasDict = BiasDictlhd(b1=1, bη=0, P0_ϵ=0, b∇2δ=mcP_Lv(Pk)[1]/6, β∇2v=0,
                                   β∂2parav=0, P2_ϵ=0, P2_ϵε_η=0)
    k, p13EFT = P13ggs(Pk, [1.0], BiasOperatorDict13=local13BiasDict, f=0)
    pL13 = 2*p13EFT .+ Plhd(Pk, k, [1.0], BiasOperatorDictlhd=locallhdBiasDict, f=0)
    p13 = (pL13[:,1] .- Pk.(k))
    return k, p13
end

function P13BiasTest2f(Pk)
    local13BiasDict = BiasDict13(b1=1, bη=0, bK²=0, bδη=0, bη²=0, bKKpara=0, btd=0,
                                 bδΠ2para=0, bηΠ2para=0, bΠ2Kpara=0, bΠ2para=0, bΠ3para=0)
    locallhdBiasDict = BiasDictlhd(b1=1, bη=0, P0_ϵ=0, b∇2δ=mcP_Lv(Pk)[1]/6, β∇2v=0,
                                   β∂2parav=0, P2_ϵ=0, P2_ϵε_η=0)
    k, p13EFT = P13ggs(Pk, [1.0], BiasOperatorDict13=local13BiasDict, f=1)
    pL13 = 2*p13EFT .+ Plhd(Pk, k, [1.0], BiasOperatorDictlhd=locallhdBiasDict, f=1)
    p13 = (pL13[:,1] .- Pk.(k))
    return k, p13
end

function Itest(Pk)
    k, P1 = I1(Pk)
    _, P2 = I2(Pk)
    _, P3 = I3(Pk)
    _, P4 = I4(Pk)
    _, P5 = I5(Pk)
    plot(k,P1, label="I1", c="r")
    plot(k,P2, label="I2", c="b", ls="-.")
    plot(k,-P3, label="-I3", c="g", ls="--")
    plot(k, P4, label="I4", c="purple")
    plot(k, -P4, label="-I4", c="purple", ls="--")
    plot(k, P5, label="I5", c="orange", ls=":")
    xscale("log")
    yscale("log")
    legend()
    xlim(1e-3,1)
    ylim(1e-7,4)
end