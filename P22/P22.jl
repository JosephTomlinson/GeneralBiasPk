function intMnli(Mnli, r, n, l, i; q="best")
    if all(Mnli[n,l,i] .== 0)
        return [0.0],[0.0]
    end
    Mfunc = Spline1D(r,Mnli[n,l,i])
    if q == "best"
        q = best_q(Mfunc, i, 0, rmini, rmaxi)
    end
    #@show q, n, l, i
    k, intM = xicalc(Mfunc, i, 0, N=N, kmin=rmini, kmax=rmaxi, r0=r0i, q=q)
    return k, intM
end

function P22lggs(Pk, l; BiasOperatorDict22, f, Mnli, r)
    P = zeros(Float64,N)
    k,_ = intMnli(Mnli, r, 0, 0, 0)
    for i = 0:8
        for n = 0:4
            try
                mint = k.^n .* intMnli(Mnli, r, n, l, i)[2]
                P .+= (mint .- mint[1])#k.^n .* intMnli(Mnli, r, n, l, i)[2]
                #@show P[1], n, l, i
            catch y
                if isa(y,UndefRefError)
                
                else
                    @show y, n, l, i
                end
            end
        end
    end
    return k, 2*(2*pi)^3 .* P
end

function P22ggsmultipole(Pk, μs; BiasOperatorDict22, f)
    L(n,x) = Legendre(n,x)
    r, Mnli = generateMnliFast(Pk, BiasOperatorDict22, f)
    Pl = l -> P22lggs(Pk, l, BiasOperatorDict22=BiasOperatorDict22, f=f, Mnli=Mnli, r=r)
    k, P0 = Pl(0)
    _, P2 = Pl(2)
    _, P4 = Pl(4)
    _, P6 = Pl(6)
    _, P8 = Pl(8)
    p22mp = Array{Float64}(undef, length(k), length(μs))
    for i = 1:length(μs)
        @. p22mp[:,i] = L(0,μs[i])*P0 + L(2,μs[i])*P2 + L(4,μs[i])*P4 + L(6,μs[i])*P6 + L(8,μs[i])*P8
    end
    return k, p22mp
end