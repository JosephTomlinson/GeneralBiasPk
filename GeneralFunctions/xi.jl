@memoize function ξ(Pk,ℓ,n; q="best") # Eq 18 in SVM + Eq 3 in 2Fast.
    N = 16*1024
    kmino = exp(-25)
    kmaxo = exp(25)  
    r0o = 1/kmaxo
    rmini = exp(-11)#exp(-14.555)
    rmaxi = exp(11)#exp(11)
    
    
    if q == "auto"
        if ℓ == 2 && n == -2
            q = 0
        else
            q = max(.5, 1+n)
        end
    end
    
    if q == "best"
        q = best_q(Pk,ℓ,n, kmino, kmaxo)
        if ℓ == 0
            q += .4
        end
        #if ℓ == 1
        #    if n == 0
        #        q -= 1
        #    end
        #end
    end
    
    if q == "show"
        q = show_best_q(Pk,ℓ,n, kmino, kmaxo)
        if ℓ == 0
            println("Xi l0 correction")
            q += .4
        end
        #if ℓ == 1
        #    if n == 0
        #        println("Xi l1n0 correction")
        #        q -= 1
        #    end
        #end
    end
    
    if q == "henry"
        q = henry_q(Pk,ℓ,n, kmino, kmaxo, r=kmaxo)#r=rmaxi)
    end
    
    r, xi = xicalc(Pk, ℓ, -n, N=N, kmin=kmino, kmax=kmaxo, r0=r0o; q=q)
    return r, xi .* (r .^ (-n))
end