function generateCln(BiasDict13, f)
    b1 = BiasDict13["δ3"]
    bη = BiasDict13["η3"]
    bK² = BiasDict13["2tr[KK2]"]
    bδη = BiasDict13["δη2"]
    bη² = BiasDict13["2ηη2"]
    bKKpara = BiasDict13["2KK2∥"]
    btd = BiasDict13["Otd"]
    bδΠ2para = BiasDict13["δΠ2∥"]
    bηΠ2para = BiasDict13["ηΠ2∥"]
    bΠ2Kpara = BiasDict13["Π2K∥"]
    bΠ2para = -BiasDict13["s^k∂kΠ2∥"]
    bΠ3para = BiasDict13["Π3∥"] - 2bΠ2para
    
    Cln = zeros(Float64, 7, 5)
    
    Cln[1,1] = (2*b1^2)/21. + (10*b1*bK²)/7. + (4*b1*btd)/7. - (8*f*b1*bη)/63. - 
    (10*f*bK²*bη)/21. - (4*f*btd*bη)/21. + (2*f^2*bη^2)/35. + (6*f^3*bη^2)/245. + 
    (4*f^2*b1*bη²)/35. - (12*f^3*bη*bη²)/245. + (10*b1*bKKpara)/21. - 
    (4*f*bη*bKKpara)/21. - (2*f*b1*bηΠ2para)/21. + (2*f^2*bη*bηΠ2para)/49. + 
    (4*b1*bΠ2para)/63. + (4*f*bη*bΠ2para)/105. - (f^2*bη*bΠ2para)/49. + 
    (5*b1*bΠ2Kpara)/21. - (2*f*bη*bΠ2Kpara)/21. + (2*b1*bΠ3para)/63. + (3*f*bη*bΠ3para)/70.
    
    Cln[1,2] = b1^2/14. + (f*b1^2)/7. - (15*b1*bK²)/7. - (6*b1*btd)/7. + (f*b1*bη)/21. - 
    (6*f^2*b1*bη)/35. + (5*f*bK²*bη)/7. + (2*f*btd*bη)/7. - (3*f^2*bη^2)/70. + 
    (6*f^3*bη^2)/245. - (6*f^2*b1*bη²)/35. + (18*f^3*bη*bη²)/245. - 
    (5*b1*bKKpara)/7. + (2*f*bη*bKKpara)/7. + (f*b1*bηΠ2para)/7. - 
    (3*f^2*bη*bηΠ2para)/49. - (2*b1*bΠ2para)/21. - (2*f*bη*bΠ2para)/35. + 
    (3*f^2*bη*bΠ2para)/98. - (5*b1*bΠ2Kpara)/14. + (f*bη*bΠ2Kpara)/7. - 
    (b1*bΠ3para)/21. - (9*f*bη*bΠ3para)/140.
    
    Cln[1,3] = -b1^2/6. - (f*b1^2)/7. + (5*b1*bK²)/7. + (2*b1*btd)/7. + (5*f*b1*bη)/63. + 
    (2*f^2*b1*bη)/35. - (5*f*bK²*bη)/21. - (2*f*btd*bη)/21. - 
    (4*f^2*bδη*bη)/35. - (f^2*bη^2)/70. + (12*f^3*bη^2)/245. + 
    (2*f^2*b1*bη²)/35. + (18*f^3*bη*bη²)/245. + (5*b1*bKKpara)/21. - 
    (2*f*bη*bKKpara)/63. + (4*f*bη*bδΠ2para)/21. - (f*b1*bηΠ2para)/21. - 
    (3*f^2*bη*bηΠ2para)/49. + (2*b1*bΠ2para)/63. + (6*f*bη*bΠ2para)/35. - 
    (5*f^2*bη*bΠ2para)/98. + (5*b1*bΠ2Kpara)/42. - (f*bη*bΠ2Kpara)/63. + 
    (b1*bΠ3para)/63. + (61*f*bη*bΠ3para)/420.
    
    Cln[1,4] = (4*f^2*b1*bη)/35. + (4*f^2*bδη*bη)/35. - (24*f^3*bη^2)/245. - 
    (24*f^3*bη*bη²)/245. - (4*f*bη*bKKpara)/63. - (4*f*bη*bδΠ2para)/21. + 
    (4*f^2*bη*bηΠ2para)/49. - (16*f*bη*bΠ2para)/105. + (2*f^2*bη*bΠ2para)/49. - 
    (2*f*bη*bΠ2Kpara)/63. - (13*f*bη*bΠ3para)/105.
    
    Cln[3,1] = (-16*f*b1*bη)/63. - (9*f^2*b1*bη)/49. - (20*f*bK²*bη)/21. - 
    (8*f*btd*bη)/21. + (8*f^2*bη^2)/49. + (f^3*bη^2)/7. + (4*f^2*b1*bη²)/49. - 
    (4*f^3*bη*bη²)/49. + (5*b1*bKKpara)/21. - (65*f*bη*bKKpara)/147. - 
    (10*f*b1*bηΠ2para)/147. + (10*f^2*bη*bηΠ2para)/147. - (4*b1*bΠ2para)/9. + 
    (15*f*b1*bΠ2para)/98. + (4*f*bη*bΠ2para)/21. - (5*f^2*bη*bΠ2para)/42. + 
    (5*b1*bΠ2Kpara)/42. - (65*f*bη*bΠ2Kpara)/294. - (101*b1*bΠ3para)/252. + 
    (37*f*bη*bΠ3para)/196.
    
    Cln[3,2] = (2*f*b1^2)/7. + (2*f*b1*bη)/21. - (3*f^2*b1*bη)/14. + (10*f*bK²*bη)/7. + 
    (4*f*btd*bη)/7. - (6*f^2*bη^2)/49. - (f^3*bη^2)/98. - (6*f^2*b1*bη²)/49. + 
    (6*f^3*bη*bη²)/49. - (5*b1*bKKpara)/14. + (65*f*bη*bKKpara)/98. + 
    (5*f*b1*bηΠ2para)/49. - (5*f^2*bη*bηΠ2para)/49. + (2*b1*bΠ2para)/3. - 
    (45*f*b1*bΠ2para)/196. - (2*f*bη*bΠ2para)/7. + (5*f^2*bη*bΠ2para)/28. - 
    (5*b1*bΠ2Kpara)/28. + (65*f*bη*bΠ2Kpara)/196. + (101*b1*bΠ3para)/168. - 
    (111*f*bη*bΠ3para)/392.
    
    Cln[3,3] = (4*f*b1^2)/7. + (6*f*b1*bδη)/7. + (10*f*b1*bη)/63. - (11*f^2*b1*bη)/14. - 
    (10*f*bK²*bη)/21. - (4*f*btd*bη)/21. - (22*f^2*bδη*bη)/49. - (2*f^2*bη^2)/49. + 
    (11*f^3*bη^2)/98. - (34*f^2*b1*bη²)/49. + (10*f^3*bη*bη²)/49. - (5*b1*bKKpara)/14. + 
    (25*f*bη*bKKpara)/882. - (10*b1*bδΠ2para)/7. + (110*f*bη*bδΠ2para)/147. + 
    (85*f*b1*bηΠ2para)/147. - (25*f^2*bη*bηΠ2para)/147. - (86*b1*bΠ2para)/63. + 
    (75*f*b1*bΠ2para)/196. + (34*f*bη*bΠ2para)/49. - (5*f^2*bη*bΠ2para)/196. - 
    (5*b1*bΠ2Kpara)/28. + (25*f*bη*bΠ2Kpara)/1764. - (569*b1*bΠ3para)/504. + 
    (683*f*bη*bΠ3para)/1176.
    
    Cln[3,4] = (-6*f*b1^2)/7. - (6*f*b1*bδη)/7. + (58*f^2*b1*bη)/49. + 
    (22*f^2*bδη*bη)/49. + (4*f^3*bη^2)/49. + (36*f^2*b1*bη²)/49. + 
    (4*f^3*bη*bη²)/49. + (10*b1*bKKpara)/21. - (110*f*bη*bKKpara)/441. + 
    (10*b1*bδΠ2para)/7. - (110*f*bη*bδΠ2para)/147. - (30*f*b1*bηΠ2para)/49. - 
    (10*f^2*bη*bηΠ2para)/147. + (8*b1*bΠ2para)/7. - (15*f*b1*bΠ2para)/49. - 
    (88*f*bη*bΠ2para)/147. - (15*f^2*bη*bΠ2para)/49. + (5*b1*bΠ2Kpara)/21. - 
    (55*f*bη*bΠ2Kpara)/441. + (13*b1*bΠ3para)/14. - (143*f*bη*bΠ3para)/294.
    
    Cln[3,5] = (-16*f^3*bη^2)/49. - (16*f^3*bη*bη²)/49. + (40*f^2*bη*bηΠ2para)/147. + 
    (40*f^2*bη*bΠ2para)/147.
    
    Cln[5,1] = (-12*f^2*b1*bη)/49. + (16*f^2*bη^2)/245. + (72*f^3*bη^2)/385. - 
    (48*f^2*b1*bη²)/245. + (192*f^3*bη*bη²)/2695. - (4*f*bη*bKKpara)/49. + 
    (8*f*b1*bηΠ2para)/49. - (32*f^2*bη*bηΠ2para)/539. + (10*f*b1*bΠ2para)/49. + 
    (16*f*bη*bΠ2para)/105. - (12*f^2*bη*bΠ2para)/77. - (2*f*bη*bΠ2Kpara)/49. + 
    (101*f*bη*bΠ3para)/735.
    
    Cln[5,2] = (6*f^2*b1*bη)/35. - (12*f^2*bη^2)/245. - (36*f^3*bη^2)/245. + 
    (72*f^2*b1*bη²)/245. - (288*f^3*bη*bη²)/2695. + (6*f*bη*bKKpara)/49. - 
    (12*f*b1*bηΠ2para)/49. + (48*f^2*bη*bηΠ2para)/539. - (15*f*b1*bΠ2para)/49. - 
    (8*f*bη*bΠ2para)/35. + (18*f^2*bη*bΠ2para)/77. + (3*f*bη*bΠ2Kpara)/49. - 
    (101*f*bη*bΠ3para)/490.
    
    Cln[5,3] = (18*f^2*b1*bη)/35. - (72*f^2*bδη*bη)/245. - (4*f^2*bη^2)/245. - 
    (432*f^3*bη^2)/2695. + (156*f^2*b1*bη²)/245. - (228*f^3*bη*bη²)/2695. + 
    (6*f*bη*bKKpara)/49. + (24*f*bη*bδΠ2para)/49. - (26*f*b1*bηΠ2para)/49. + 
    (38*f^2*bη*bηΠ2para)/539. - (45*f*b1*bΠ2para)/49. + (344*f*bη*bΠ2para)/735. + 
    (180*f^2*bη*bΠ2para)/539. + (3*f*bη*bΠ2Kpara)/49. + (569*f*bη*bΠ3para)/1470.
    
    Cln[5,4] = (-528*f^2*b1*bη)/245. + (72*f^2*bδη*bη)/245. + (2664*f^3*bη^2)/2695. - 
    (120*f^2*b1*bη²)/49. + (2664*f^3*bη*bη²)/2695. - (8*f*bη*bKKpara)/49. - 
    (24*f*bη*bδΠ2para)/49. + (100*f*b1*bηΠ2para)/49. - (444*f^2*bη*bηΠ2para)/539. + 
    (120*f*b1*bΠ2para)/49. - (96*f*bη*bΠ2para)/245. - (612*f^2*bη*bΠ2para)/539. - 
    (4*f*bη*bΠ2Kpara)/49. - (78*f*bη*bΠ3para)/245.
    
    Cln[5,5] = (12*f^2*b1*bη)/7. - (468*f^3*bη^2)/539. + (12*f^2*b1*bη²)/7. - 
    (468*f^3*bη*bη²)/539. - (10*f*b1*bηΠ2para)/7. + (390*f^2*bη*bηΠ2para)/539. - 
    (10*f*b1*bΠ2para)/7. + (390*f^2*bη*bΠ2para)/539.
    
    Cln[7,1] = (40*f^3*bη^2)/539. + (32*f^3*bη*bη²)/539. - (80*f^2*bη*bηΠ2para)/1617. - 
    (100*f^2*bη*bΠ2para)/1617.
    
    Cln[7,2] = (-4*f^3*bη^2)/49. - (48*f^3*bη*bη²)/539. + (40*f^2*bη*bηΠ2para)/539. + 
    (50*f^2*bη*bΠ2para)/539.
    
    Cln[7,3] = (-116*f^3*bη^2)/539. - (104*f^3*bη*bη²)/539. + (260*f^2*bη*bηΠ2para)/1617. + 
    (150*f^2*bη*bΠ2para)/539.
    
    Cln[7,4] = (400*f^3*bη^2)/539. + (400*f^3*bη*bη²)/539. - (1000*f^2*bη*bηΠ2para)/1617. - 
    (400*f^2*bη*bΠ2para)/539.
    
    Cln[7,5] = (-40*f^3*bη^2)/77. - (40*f^3*bη*bη²)/77. + (100*f^2*bη*bηΠ2para)/231. + 
    (100*f^2*bη*bΠ2para)/231.
    
    return Cln
end