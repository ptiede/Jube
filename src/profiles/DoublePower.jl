"""
    DoublePower
A double power law for a profile
    ``\\frac{r}{r_0}^a\\left(1+\\frac{r}{r_0}^{-(a+b)}\\right)^{-1}``
"""
struct DoublePower <: AbstractProfile
    """ FluidVelocity struct"""
    fluid_velocity::FluidVelocity{Float64}
    """ MagneticField struct"""
    magnetic_field::MagneticField{Float64}
    """emission spectral index"""
    α::Float64
    """emission cross product spectral index"""
    αζ::Float64
    """Function that computes the profile as radius `r` and redshift `z`"""
    r0::Float64
    p1::Float64
    p2::Float64
end
@inline profile(p::DoublePower, r) = return (r/p.r0)^p.p1/(1 + (r/p.r0)^(p.p1+p.p2))

