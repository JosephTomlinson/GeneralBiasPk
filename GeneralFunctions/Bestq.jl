using Combinatorics
function henry_q(Pk,l,n,kmin,kmax; r)
    function I_m(m, smin, smax, kmin, kmax, pmin, pmax)
        kturn = (kmin^smin * kmax^(-smax) * abs(pmax/pmin))^(1/(smin-smax))
        lo = kmin^(-smin) * pmin * (kturn^(m+smin) - kmin^(m+smin))/(m+smin)
        hi = kmax^(-smax) * pmax * (kmax^(m+smax) - kmin^(m+smax))/(m+smax)
        return lo + hi
    end
    pmin = Pk(kmin)
    pmax = Pk(kmax)
    
    kvals = exp10.(range(log10(kmin), stop=log10(kmax), length=1024))
    if Pk(kvals[1]) > 0
        smin = (log(Pk(kvals[20])) - log(Pk(kvals[1]))) / (log(kvals[20]) - log(kvals[1]))
    else
        smin = (log(abs(Pk(kvals[20]))) - log(abs(Pk(kvals[1])))) / (log(kvals[20]) - log(kvals[1]))
    end
    if Pk(kvals[end]) > 0
        smax = (log(Pk(kvals[end])) - log(Pk(kvals[end-19]))) / (log(kvals[end]) - log(kvals[end-19]))
    else
        smax = (log(abs(Pk(kvals[end]))) - log(abs(Pk(kvals[end-19])))) / (log(kvals[end]) - log(kvals[end-19]))
    end
    
    Ia = I_m(3+l+n, smin, smax, kmin, kmax, pmin, pmax)
    Ib = I_m(1+n, smin, smax, kmin, kmax, pmin, pmax)
    G = log(kmax/kmin)
    qnu = 1 - l/2 + (1/2G)*log((r^(l+2)/Int64(doublefactorial(2l+1)))*(Ia/Ib))
    #@show (r^(l+1)/Int64(doublefactorial(2l+1)))*(Ia/Ib)
    qmin = max(smax+4-1+n,-l+.5)
    qmax = 3+smin+n #min(3+n1+n,2)
    
    if qnu > qmax
        q = qmax
    elseif qnu < qmin
        q = qmin
    else
        q= qnu
    end
    
    return q#+.5
end
    

function best_q(Pk,ℓ,n, kmin, kmax)
    kvals = exp10.(range(log10(kmin), stop=log10(kmax), length=1024))
    #if Pk(kvals[1]) > 0
    #    smin = (log(Pk(kvals[20])) - log(Pk(kvals[1]))) / (log(kvals[20]) - log(kvals[1]))
    #else
        smin = (log(abs(Pk(kvals[20]))) - log(abs(Pk(kvals[1])))) / (log(kvals[20]) - log(kvals[1]))
    #end
    #if Pk(kvals[end]) > 0
    #    smax = (log(Pk(kvals[end])) - log(Pk(kvals[end-19]))) / (log(kvals[end]) - log(kvals[end-19]))
    #else
        smax = (log(abs(Pk(kvals[end]))) - log(abs(Pk(kvals[end-19])))) / (log(kvals[end]) - log(kvals[end-19]))
    #end
    
    #@show smin
    #@show smax
    n1 = smin
    n2 = smax+4
    qmin = max(n2-1+n,-ℓ+.5)
    qmax = 3+n1+n#min(3+n1+n,2)
    qbest = n - (smin+smax)/2
    
    #@show qbest
    #@show qmin
    #@show qmax
    if qbest > qmax
        q = qmax
    elseif qbest < qmin
        q = qmin
    else
        q=qbest
    end
    
    if sign(Pk(kvals[1])) != sign(Pk(kvals[end]))
        #println("Changing Sign correction")
        if Pk(kvals[1]) > 0
            q -= .5
        else
            q += .5
        end
    end
    #if (smin < -1e-10) && (smax < 0)
    if (smin < -1e-4) && (smax < 0)
        #println("Montonic Decreasing correction")
        if smin > -.3
            #println("Almost Flat correction")
            q -= 1.1
        else
            #q -= .2
            q -= .4
        end
    end
    
    if smax < -6.5
        #println("Rapid Decrease correction")
        q -= 1
    #elseif (smax-smin) < -5.5
    #    println("Large Second Derivative correction")
    #    q += .4
    end
    
    if (smin > 3.5) && (smax < -3.5)
        #println("Both Slopes Large correction")
        q += 3
    end
    
    if (smin < .5) && (smax < -3.5) && (smin > 0)
        #println("Flat to Quick Fall correction")
        q -= 1.2
    end
    
    #@show q
    return q
end

