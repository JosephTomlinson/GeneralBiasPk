function M13(biasop::String)
    if biasop == "Otd"
        m = [4 -6 2 0 0;
             0  0 0 0 0;
             0  0 0 0 0]
    elseif biasop == "δΠ2∥"
        m = [0 0  5  -5 0;
             0 0 -15 15 0;
             0 0  0  0  0]
    elseif biasop == "ηΠ2∥"
        m = [ 0   0   -15/4  15/2  -15/4;
             -5  15/2   20   -60    75/2;
              5 -15/2 -65/4 125/2 -175/4]
    elseif biasop == "Π2K∥"
        m = [5/4 -15/8 35/24 -5/6 0;
             5/4 -15/8 -15/8 5/2  0;
              0     0     0   0   0]
    elseif biasop == "u2∥∂∥δ"
        m = [0  0  3 -3 0;
             0 -3 -6  9 0;
             0  0  0  0 0]
    elseif biasop == "u2∥∂∥η"
        m = [0      0   -9/4   9/2  -9/4;
            -9/4  27/8  63/8 -63/2  45/2;
            15/4 -21/8 -39/8   30 -105/4]
    elseif biasop == "s^k∂kΠ2∥"
        m = [  5/4 -15/8  25/8 -5/2 0;
             -15/4  45/8 -75/8 15/2 0;
                0     0     0    0  0]
    elseif biasop == "u∥∂∥Π2∥"
        m = [   0     0    15/4 -15/2  15/4;
              15/4 -45/8 -225/8 135/2 -75/2;
             -25/4  75/8  225/8  -75  175/4]
    elseif biasop == "Π3∥"
        m = [  13/8  -39/16   65/16 -13/4 0;
             -101/24 101/16 -569/48  39/4 0;
                 0      0       0      0  0]
    elseif biasop == "δ3"
        m = [2/3 1/2 -7/6 0 0;
              0   0    0  0 0;
              0   0    0  0 0]
    elseif biasop == "η3"
        m = [0  0    0  0 0;
             -2 3/2 1/2 0 0;
             0  0    0  0 0]          
    elseif biasop == "2tr[KK2]"
        return (5/2) .* M13("Otd")
    elseif biasop == "δη2"
        return -(3/5) .* M13("δΠ2∥")
    elseif biasop == "u∥∂∥η2"
        return -(3/5) .* M13("u∥∂∥Π2∥")
    elseif biasop == "2ηη2"
        return -(6/5) .* M13("ηΠ2∥")
    elseif biasop == "2KK2∥"
        return 2 .* M13("Π2K∥")
    else
        #println("Error $biasop")
    end
    return m ./ 7
end