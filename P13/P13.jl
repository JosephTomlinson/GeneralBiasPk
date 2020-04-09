GrowthPowerDict13 = Dict([("δ3", 0), ("η3", 1) , ("2tr[KK2]", 0), ("δη2", 1),
                           ("2ηη2", 2), ("2KK2∥", 0), ("Otd", 0), ("δΠ2∥", 0),
                           ("ηΠ2∥", 1), ("Π2K∥", 0), ("s^k∂kΠ2∥", 0),
                           ("u2∥∂∥δ", 1), ("u2∥∂∥η", 2), ("u∥∂∥η2", 2), 
                           ("u∥∂∥Π2∥", 1), ("Π3∥", 0)])

@memoize function P_13ln(Pk,ℓ,n)
    r, ξℓn = ξ(Pk,ℓ,n)
    
    function get_int_func(f)
        Spline1D(r, f./r, bc="extrapolate")
    end
    if n == 0
        q = .5
    elseif n == 1
        q = .1
    elseif n == -1
        q = .5
    elseif n == -2
        q = 1.1
    elseif (ℓ == 2) && (n == 2)
        q = -.8
    elseif n == 2
        q = -1
    end
    k, mcP = xicalc(get_int_func(ξℓn), ℓ, 0, kmin=rmini, kmax=rmaxi, r0=r0i, N=N, q=q)
    return k, 2π^2 .* mcP
end

@memoize function I1(Pk)
    k, mcP1m1 = P_13ln(Pk,1,-1)
    _, mcP3m1 = P_13ln(Pk,3,-1)
    mcI1 = @. (2k^3/5)*(mcP1m1 - mcP3m1)
    return k, mcI1
end

@memoize function I2(Pk)
    k, mcP00 = P_13ln(Pk,0,0)
    _, mcP20 = P_13ln(Pk,2,0)
    mcI2 = @. (2k^2/3)*(mcP00 - mcP20)
    return k, mcI2
end

@memoize function I3(Pk)
    k, mcP02 = P_13ln(Pk,0,2)
    _, mcP22 = P_13ln(Pk,2,2)
    mcI3 = @. (2/3)*(mcP02 - mcP22)
    return k, mcI3
end

@memoize function I4(Pk)
    k, mcP02 = P_13ln(Pk,0,2)
    _, mcP22 = P_13ln(Pk,2,2)
    _, mcP42 = P_13ln(Pk,4,2)
    mcI4 = @. (2/15)*mcP02 + (2/21)*mcP22 - (8/35)*mcP42
    return k, mcI4
end

@memoize function I5(Pk)
    k, mcP02 = P_13ln(Pk,0,2)
    _, mcP22 = P_13ln(Pk,2,2)
    _, mcP42 = P_13ln(Pk,4,2)
    _, mcP62 = P_13ln(Pk,6,2)
    mcI5 = @. (2/35)*mcP02 + (2/21)*mcP22 - (32/385)*mcP42 - (16/231)*mcP62
    return k, mcI5
end

function fO(Pk,O::String,μs)
    m = M13(O)
    k, i1 = I1(Pk)
    _, i2 = I2(Pk)
    _, i3 = I3(Pk)
    _, i4 = I4(Pk)
    _, i5 = I5(Pk)
    fout = Array{Float64}(undef,length(k),length(μs))
    for j = 1:length(μs)
        for i = 1:length(k)
            fout[i,j] = ([1 μs[j]^2 μs[j]^4] * m * [i1[i], i2[i], i3[i], i4[i], i5[i]])[1]
        end
    end
    return k, fout
end

function P13ggs(Pk, μs; BiasOperatorDict13, f)
    p13 = zeros(Float64, N, length(μs))
    kout = Array{Float64}(undef, N)
    i = 1
    b1 = BiasOperatorDict13["δ3"]
    bη = BiasOperatorDict13["η3"]
    local k, fo
    for op in keys(BiasOperatorDict13)
        try
            k, fo = fO(Pk,op,μs)
        catch
            continue
        end
        
        if i == 1
            kout .= k
            i += 1
        end
        bo = BiasOperatorDict13[op]
        nf = GrowthPowerDict13[op]
        for j = 1:length(μs)
            p13[:,j] .+= (b1-bη*f*μs[j]^2) * bo * f^nf .* fo[:,j] .* Pk.(k)
        end
    end
    return kout, p13
end

function P13lggs(Pk, l; BiasOperatorDict13, f, Cln)
    if l == 8
        return [0.0], [0.0]
    end
    k, i1 = I1(Pk)
    _, i2 = I2(Pk)
    _, i3 = I3(Pk)
    _, i4 = I4(Pk)
    _, i5 = I5(Pk)
    In = [i1,i2,i3,i4,i5]
    Pkvals = Pk.(k)
    p13l = zeros(Float64, length(k))
    for n = 1:5
        p13l .+= Cln[l+1,n] .* In[n] .* Pkvals
    end
    return k, p13l
end

function P13ggsmultipole(Pk, μs; BiasOperatorDict13, f)
    L(n,x) = Legendre(n,x)
    Cln = generateCln(BiasOperatorDict13, f)
    Pl = l -> P13lggs(Pk,l, BiasOperatorDict13=BiasOperatorDict13, f=f, Cln=Cln)
    k, P0 = Pl(0)
    _, P2 = Pl(2)
    _, P4 = Pl(4)
    _, P6 = Pl(6)
    p13mp = Array{Float64}(undef, length(k), length(μs))
    for i = 1:length(μs)
        @. p13mp[:,i] = L(0,μs[i])*P0 + L(2,μs[i])*P2 + L(4,μs[i])*P4 + L(6,μs[i])*P6
    end
    return k, p13mp
end