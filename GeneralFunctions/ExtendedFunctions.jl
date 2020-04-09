struct extended_function
    fmin::Function
    spline::Function
    fmax::Function
    xmin::Float64
    xmax::Float64
end

function extend_func(k, Pk, kmin, kmax) # Needs serious modification, Log/Lin??, slope eval != kmin/kmax
    filt = (k .>= kmin) .&  (k .<= kmax) .& (Pk .> 0)
    spl = Spline1D(log.(k[filt]), log.(Pk[filt]), k=2)
    splf = x -> exp(spl(log(x)))
    
    smin = derivative(spl, log(kmin))
    smax = derivative(spl, log(kmax))
    
    fmin = x -> splf(kmin) * (x/kmin)^smin
    fmax = x -> splf(kmax) * (x/kmax)^smax

    return extended_function(fmin, splf, fmax, kmin, kmax)
end

function eval_extend_func(Pkfunc, x)
    if x < Pkfunc.xmin
        Pkfunc.fmin(x)
    elseif x < Pkfunc.xmax
        Pkfunc.spline(x)
    else
        Pkfunc.fmax(x)
    end
end

(pk::extended_function)(x) = eval_extend_func(pk,x)