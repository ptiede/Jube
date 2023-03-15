struct JuKeBOX <: AbstractEmissionModel
    nmax::Int
    metric::Union{Kerr, Schwarzschild}
    profile::DoublePower
    function JuKeBOX(nmax, spin, α, αζ, rpeak, p1, p2, βv, χ, ι, η)
        met = Kerr(spin)
        fluid_vel = FluidVelocity(βv, χ)
        mag_field = MagneticField(ι, η)
        prof = DoublePower(fluid_vel, mag_field, α, αζ, rpeak, p1, p2)
        return new(nmax, met, prof)
    end
    function JuKeBOX(nmax, α, αζ, rpeak, p1, p2, βv, χ, ι, η)
        met = Schwarzschild()
        fluid_vel = FluidVelocity(βv, χ)
        mag_field = MagneticField(ι, η)
        prof = DoublePower(fluid_vel, mag_field, α, αζ, rpeak, p1, p2)
        return new(nmax, met, prof)
    end

end

function raytrace(acc::JuKeBOX, α, β, θs, o::AssymptoticObserver, isindir::Bool)
    stokesi = 0.; stokesq = 0.; stokesu = 0.
    for n in 0:acc.nmax
        _, i, q, u = trace_nring(acc, n, α, β, θs, o, isindir)
        stokesi += i
        stokesq += q
        stokesu += u
    end

    return @SVector[stokesi, stokesq, stokesu, 0.0]#zero(typeof(stokesi)))
end

function raytrace_n(acc::JuKeBOX, n, α, β, θs, o::AssymptoticObserver, isindir)
    _, stokesi, stokesq, stokesu = trace_nring(acc, n, α, β, θs, o, isindir)

    return @SVector[stokesi, stokesq, stokesu, zero(typeof(stokesi))]
end

function raytrace_and_get_mask(acc::JuKeBOX ,n, α, β, θs, o::AssymptoticObserver, isindir)
   (_, stokesi, stokesq, stokesu), mask = trace_nring_and_get_mask(acc, n, α, β, θs, o, isindir)

    return @SVector[stokesi, stokesq, stokesu, zero(typeof(stokesi))], mask
end

function trace_nring(acc::JuKeBOX, n::Int, α, β, θs, o::AssymptoticObserver, isindir)
    met = acc.metric
    rh = horizon(met)
    r, νr, _ = rs(met, α, β, θs, o, isindir, n) 
    T = eltype(r)

    νθ =  cos(θs)< abs(cos(o.inclination)) ? (o.inclination>θs) ⊻ (n%2==1) : !isindir

    # bail out if r is within horizon
    (r == Inf || r < rh + eps()) && return zero(T), zero(T), zero(T), zero(T)

    return _emission(acc, α, β, r, νr, νθ, o, θs)
end

function trace_nring_and_get_mask(acc::JuKeBOX, n, α, β, θs, o::AssymptoticObserver, isindir)
    met = acc.metric
    rh = horizon(met)
    (r, νr, _), mask = rs_mask(met, n, α, β, θs, o, isindir) 
    T = eltype(r)

    νθ =  cos(θs)< abs(cos(o.inclination)) ? (o.inclination>θs) ⊻ (n%2==1) : !isindir

    (r == Inf || r < rh + eps()) && return (zero(T), zero(T), zero(T), zero(T)), mask

    return _emission(acc, α, β, r, νr, νθ, o, θs), mask
end

function _emission(acc::JuKeBOX, α, β, r, νr, νθ, o, θs)
    met = acc.metric
    fluid_vel = acc.profile.fluid_velocity
    fluidβ = @SVector [fluid_vel.β, π/2. , fluid_vel.χ]
    mag_field = acc.profile.magnetic_field
    mag_vec = @SVector [mag_field.br, mag_field.bϕ, mag_field.bz]
    eα, eβ, redshift, lp = calcPol(met, α, β, r, θs, o.inclination, acc.profile.αζ, mag_vec, fluidβ , νr, νθ)

    # Get the profile value at the emission radius
    prof = profile(acc.profile, r)*max(redshift , eps())^(3+acc.profile.α)
    # We add a small perturbation to q and u. This prevent taking a
    # derivative of hypot at (0,0) which gives  a NaN. This is silly because this
    # occurs when there is no emission so the derivative should be zero there.
    q = -(eα^2 - eβ^2)*lp*prof + eps()
    u = -2*eα*eβ*lp*prof + eps()
    i = hypot(q, u)
    return redshift, i, q, u
end

@inline function get_nturns(n, sβ)
    if sβ > 0
        return sβ+n
    else
        return sβ+n+1
    end
end