function Legendre(n,x)
    if n == 0
        1
    elseif n == 1
        x
    elseif n == 2
        (1/2)*(3x^2-1)
    elseif n == 3
        (1/2)*(5x^3-3x)
    elseif n == 4
        (1/8)*(35x^4-30x^2 + 3)
    elseif n == 5
        (1/8)*(63x^5 - 70x^3 + 15)
    elseif n == 6
        (1/16)*(231x^6 - 315x^4 + 105x^2 - 5)
    elseif n == 7
        (1/16)*(-35x + 315x^3 - 693x^5 + 429x^7)
    elseif n == 8
        (1/128)*(35 - 1260x^2 + 6930x^4 - 12012x^6 + 6435x^8)
    end
end