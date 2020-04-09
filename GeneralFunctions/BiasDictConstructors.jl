function BiasDict13(; b1, bη, bK², bδη, bη², bKKpara, btd, bδΠ2para, bηΠ2para, bΠ2Kpara, bΠ2para, bΠ3para)
    BiasOperatorDict13 = Dict([("δ3", b1), ("η3", bη) , ("2tr[KK2]", bK²), ("δη2",bδη),
                           ("2ηη2", bη²), ("2KK2∥", bKKpara), ("Otd", btd), ("δΠ2∥", bδΠ2para),
                           ("ηΠ2∥", bηΠ2para), ("Π2K∥", bΠ2Kpara), ("s^k∂kΠ2∥", -bΠ2para),
                           ("u2∥∂∥δ", -b1), ("u2∥∂∥η", -bη), ("u∥∂∥η2", -bη), 
                           ("u∥∂∥Π2∥", -bΠ2para), ("Π3∥", bΠ3para+2bΠ2para)])
end

function BiasDictlhd(; b1, bη, P0_ϵ, b∇2δ, β∇2v, β∂2parav, P2_ϵ, P2_ϵε_η)
    BiasOperatorDictlhd = Dict([("b1", b1),("bη", bη), ("P0ϵ", P0_ϵ), ("b∇2δ", b∇2δ), ("β∇2v", β∇2v), 
                            ("β∂2∥v", β∂2parav), ("P2ϵ", P2_ϵ), ("P2ϵεη", P2_ϵε_η)])
end

function BiasDict22(; b1, bη, b2, bK², bδη, bη², bKKpara, bΠ2para)
    BiasOperatorDict22 = Dict([("b1",b1), ("bη", bη), ("b2", b2), ("bK2", bK²), ("bδη", bδη),
                           ("bη2", bη²), ("bKK∥", bKKpara), ("bΠ2∥", bΠ2para)])
end