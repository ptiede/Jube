"""
    GaussianRing
A gaussian ring profile
"""
struct GaussianRing <: AbstractProfile
    """ FluidVelocity struct"""
    fluid_velocity::FluidVelocity{Float64}
    """ MagneticField struct"""
    magnetic_field::MagneticField{Float64}
    """emission spectral index"""
    α::Float64
    """emission cross product spectral index"""
    αζ::Float64
    """Function that computes the profile as radius `r` and redshift `z`"""
    rpeak::Float64
    width::Float64
end
@inline profile(p::GaussianRing, r) = exp(-4.0*log(2.0)*abs2( (r-p.rpeak)/p.width ))

