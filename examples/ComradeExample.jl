import Pkg;
Pkg.activate((@__DIR__))
using Comrade
using Plots
using Jube

params  = NamedTuple{(:m_d, :spin, :θo, :pa, :spec_index, :cross_spec_index, :rpeak, :p1, :p2, :βv, :χ, :ι)}([3.7855118062030972, -0.8315634748162867, 9.419251185101801*π/180, 3.161924161899382 , 0.23968279220600633, -1.4554064047766926, 4.685817341405595, 1.4050257636562866, 5.69456394172345, 0.973957150235823, -2.187948907269933, 1.3113050873728065])
nmax = 1

abstract type AccretionModel <: ComradeBase.AbstractModel end
ComradeBase.visanalytic(::Type{<:AccretionModel}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:AccretionModel}) = ComradeBase.IsAnalytic()

struct JKConeModel{B,S,O} <: AccretionModel
    acc::B
    s::S
    o::O
end

ComradeBase.isprimitive(::Type{<:JKConeModel}) = ComradeBase.IsPrimitive()

function ComradeBase.intensity_point(s::JKConeModel, p)
    α=p.X
    β=p.Y
    direct  = Jube.raytrace(s.acc, -α, β, s.s, s.o, true)[1]
    indirect = Jube.raytrace(s.acc, -α, β, s.s, s.o, false)[1]

    return direct + indirect
    
end

struct SModel{C}
    cache::C
end
function(Smodel::SModel)(θ)
    (;m_d, spin, θo, rpeak, pa, p1, p2, χ, ι, βv, spec_index, cross_spec_index) = θ
    θs = π/2
    acc = Jube.JuKeBOX(nmax, spin, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
    observer = Jube.AssymptoticObserver(1, θo)
    m = JKConeModel(acc, θs, observer)
    mm = modify(m, Stretch(μas2rad(m_d), μas2rad(m_d)), Rotate(pa))
    mimg = Comrade.modelimage(mm, Smodel.cache)
    return mimg
end


sze = 120
modelfov = 80
Smodel = SModel(create_cache(FFTAlg(padfac=8), IntensityMap(zeros(sze, sze), μas2rad(modelfov), μas2rad(modelfov))))

smodel = Smodel(params)
Plots.plot(Smodel(params))