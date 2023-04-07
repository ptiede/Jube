module Jube
import FastElliptic
using DocStringExtensions
using StaticArrays
import Elliptic
using ArbNumerics
using Plots
# Profile Exports
#--------------------------------------------------------------------------
struct FluidVelocity{T}
    β::T
    χ::T
    βr::T
    βϕ::T
    function FluidVelocity(β, χ)
        T = promote_type(eltype(β), eltype(χ))
        βT, χT = promote(β, χ)
        s,c = sincos(χT)
        return new{T}(βT, χT, βT*c, βT*s)
    end
end
struct MagneticField{T}
    ι::T
    η::T
    beq::T
    bz::T
    br::T
    bϕ::T
    function MagneticField(ι, η)
        T = promote_type(typeof(ι), typeof(η))
        ιT, ηT = promote(ι, η)
        beq,bz = sincos(ι)
        s,c = sincos(η)
        br = beq*c
        bϕ = beq*s
        return new{T}(ιT, ηT, beq, bz, br, bϕ)
    end
end
abstract type AbstractProfile end
function profile(p::AbstractProfile, args...) @error("Profile is not defined for $(typeof(p))") end

export DoublePower
include("profiles/DoublePower.jl")
export GaussianRing
include("profiles/GaussianRing.jl")

# Metric Exports
#--------------------------------------------------------------------------
export metric_uu, metric_dd, horizon, Kerr, get_roots, rs, θs, calcPol, η, λ

# Space Time Event Declaration
export AssympototicObserver, get_coordinates

abstract type AbstractSpaceTimeEvent end
function get_coordinates(event::AbstractSpaceTimeEvent) 
    throw("get_coordinates has not been defined for even of type: $(typeof(event))")
end
struct AssymptoticObserver{P,O} <: AbstractSpaceTimeEvent
    azimuth::P
    inclination::O
end
function get_coordinates(observer::AssymptoticObserver)
    return (azimuth=observer.azimuth, inclination=observer.inclination)
end

abstract type AbstractMetric end
#function metric_dd(met::AbstractMetric) @error("The metric has not been defined for $(typeof(met))")end
function metric_dd(metric::AbstractMetric, args...) ::AbstractMatrix
    throw(MethodError(metric_dd, (metric, args...)))
    return [1 0; 0 1]
end
function metric_uu(metric::AbstractMetric, args...) ::AbstractMatrix
    try return inv(metric_dd(metric::AbstractMetric, args...)) catch err throw(err) end 
end
include("metrics/Kerr.jl")

# Emission Exports
#--------------------------------------------------------------------------
abstract type AbstractEmissionModel end
function raytrace(acc::AbstractEmissionModel, varargs...)
    @error "raytrace has not been defined for raytrace($(typeof(acc)))"
end

export JuKeBOX
include("emissionModels/JuKeBOX.jl")

end
