import Pkg;
Pkg.activate((@__DIR__))
using Comrade
using Plots
using Chain
import JuKeBOX as JK
include((@__DIR__)*"/src/Jube.jl")
include((@__DIR__)*"/JKCone.jl")

params  = NamedTuple{(:m_d, :spin, :θo, :pa, :spec_index, :cross_spec_index, :rpeak, :p1, :p2, :βv, :χ, :ι)}([3.7855118062030972, -0.8315634748162867, 9.419251185101801*π/180, 3.161924161899382 , 0.23968279220600633, -1.4554064047766926, 4.685817341405595, 1.4050257636562866, 5.69456394172345, 0.973957150235823, -2.187948907269933, 1.3113050873728065])
nmax = 1

abstract type AccretionModel <: ComradeBase.AbstractPolarizedModel end
ComradeBase.visanalytic(::Type{<:AccretionModel}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:AccretionModel}) = ComradeBase.IsAnalytic()

struct ConeModel{B,S,O} <: AccretionModel
    acc::B
    s::S
    o::O
end

ComradeBase.isprimitive(::Type{<:ConeModel}) = ComradeBase.IsPrimitive()

function ComradeBase.intensity_point(s::ConeModel, α, β)
    direct  = Jube.raytrace(s.acc, -α, β, s.s, s.o, true)
    indirect = Jube.raytrace(s.acc, -α, β, s.s, s.o, false)
    #intensity value and masks
    return direct[1] + indirect[1]#, direct[2][2]+indirect[2][2]
    
end

struct SModel{C}
    cache::C
end
function(Smodel::SModel)(θ)
    (;m_d, θo, rpeak, pa, p1, p2, χ, ι, βv, spec_index, cross_spec_index) = θ
    θs = π/2
    #pa = 0.0
    acc = Jube.JuKeBOX(nmax, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
    observer = Jube.AssymptoticObserver(1, θo)
    m = ConeModel(acc, θs, observer)
    mm = @chain m begin
        stretched(_, μas2rad(m_d), μas2rad(m_d))
        rotated(_, -pa)
    end
    mimg = Comrade.modelimage(mm, Smodel.cache, ThreadedEx())
    return mimg
end


struct Model{C}
    cache::C
end
function(model::Model)(θ)
    (;m_d, θo, rpeak, spin, pa, p1, p2, χ, ι, βv, spec_index, cross_spec_index) = θ
    θs = π/2
    #pa = 0.0
    acc = Jube.JuKeBOX(nmax, spin, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
    #g, acc = bamDblPower(nmax, spin, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
    observer = Jube.AssymptoticObserver(1, θo)
    m = ConeModel(acc, θs, observer)
    mm = @chain m begin
        stretched(_, μas2rad(m_d), μas2rad(m_d))
        rotated(_, -pa)
    end
    mimg = Comrade.modelimage(mm, model.cache, ThreadedEx())
    return mimg
end

struct JKModel{C}
    cache::C
end

function(JKmodel::JKModel)(θ)
    (;m_d, θo, rpeak, spin, pa, p1, p2, χ, ι, βv, spec_index, cross_spec_index) = θ
    θs = π/2
    #pa = 0.0
    g, acc = JK.bamDblPower(nmax, spin, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
    m = JKConeModel(acc, g, θs, θo)

    mm = @chain m begin
        stretched(_, μas2rad(m_d), μas2rad(m_d))
        rotated(_, -pa)
    end
    return Comrade.modelimage(mm, JKmodel.cache, ThreadedEx())
end

sze = 120
modelfov = 80
model = Model(create_cache(FFTAlg(padfac=8), IntensityMap(zeros(sze, sze), μas2rad(modelfov), μas2rad(modelfov))))
JKmodel = JKModel(create_cache(FFTAlg(padfac=8), IntensityMap(zeros(sze, sze), μas2rad(modelfov), μas2rad(modelfov))))
Smodel = SModel(create_cache(FFTAlg(padfac=8), IntensityMap(zeros(sze, sze), μas2rad(modelfov), μas2rad(modelfov))))

Smodel(params)
Plots.plot(model(params))
Plots.plot(JKmodel(params))
Plots.plot(Smodel(params))
model_image = model(params)
JKmodel_image = JKmodel(params)
Plots.plot(intensitymap(model_image, μas2rad(modelfov), μas2rad(modelfov), sze, sze) - intensitymap(JKmodel_image, μas2rad(modelfov), μas2rad(modelfov), sze, sze))

import CSV
using DataFrames
df = CSV.File("/Volumes/Working/Documents/Harvard/Physics/Research/EHT/bamextension/test/grmhd_bam_jk_vis_rotated/JuKeBOX_dynesty_chain.csv") |> DataFrame
test_df = df[begin:100:end, :]
max_dif = 0
#for test_case in eachrow(test_df)
#    params = test_case |> NamedTuple
#
#    tempmax = maximum(intensitymap(model(params), μas2rad(modelfov), μas2rad(modelfov), sze, sze) - intensitymap(JKmodel(params), μas2rad(modelfov), μas2rad(modelfov), sze, sze)) 
#    maxdif = tempmax > max_dif ? tempmax : max_dif
#end
#
#max_dif
#using BenchmarkTools
#@btime intensitymap(model_image, μas2rad(modelfov), μas2rad(modelfov), sze, sze)
#@btime intensitymap(JKmodel_image, μas2rad(modelfov), μas2rad(modelfov), sze, sze)