import JuKeBOX as JK

abstract type AccretionModel <: ComradeBase.AbstractPolarizedModel end
ComradeBase.visanalytic(::Type{<:AccretionModel}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:AccretionModel}) = ComradeBase.IsAnalytic()

struct JKConeModel{B,G,S,O} <: AccretionModel
    acc::B
    g::G
    s::S
    o::O
end

ComradeBase.isprimitive(::Type{<:JKConeModel}) = ComradeBase.IsPrimitive()

function ComradeBase.intensity_point(s::JKConeModel, α, β)
    direct  = JK.raytrace(-α, β, s.g, s.s, JK.Observer(1., s.o), s.acc, true)
    indirect = JK.raytrace(-α, β, s.g, s.s, JK.Observer(1., s.o), s.acc, false)
    #intensity value and masks
    return direct[1] + indirect[1]#, direct[2][2]+indirect[2][2]
    
end
