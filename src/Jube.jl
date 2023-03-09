module Jube
import FastElliptic
using DocStringExtensions
using StaticArrays

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
# Space Time Event Declaration
abstract type AbstractSpaceTimeEvent end
function getCoordinates(event::AbstractSpaceTimeEvent) 
    @error "getCoordinates has not been defined for $(typeof(event))"
end

struct AssymptoticObserver{P,O} <: AbstractSpaceTimeEvent
    azimuth::P
    inclination::O
end
function getCoordinates(observer::AssymptoticObserver)
    return (distance=observer.distance, observer=observer.inclination)
end

abstract type AbstractMetric end
function met_dd(met::AbstractMetric, args...) @error("The metric has not been defined for $(typeof(met))")end
function met_uu(met::AbstractMetric, args...) inv(met_dd(met::AbstractMetric, args...)) end
#TODO: Define a default inversion scheme

export Kerr
include("metrics/Kerr.jl")
export Schwarzschild
include("metrics/Schwarzschild.jl")
export Minkowski
include("metrics/Minkowski.jl")

# Emission Exports
#--------------------------------------------------------------------------
abstract type AbstractEmissionModel end
function raytrace(acc::AbstractEmissionModel, varargs...)
    @error "raytrace has not been defined for raytrace($(typeof(acc)))"
end

export JuKeBOX
include("emissionModels/JuKeBOX.jl")

end
