# Metric Definition
##----------------------------------------------------------------------------------------------------------------------
struct Minkowski <: AbstractMetric end

"""
    met_uu(r, θ, a)

Inverse Minkowski Metric in Boyer Lindquist (BL) coordinates.

    `met` : Minkowski Metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function met_uu(met::Minkowski, r, θ)
  r_inv = 1/r
  r2_inv = r_inv^2

  return @SMatrix [
    1   0.0   0.0     0.0
    0.0   1.0     0.0     0.0
    0.0   0.0   r2_inv  0.0
    0.0   0.0   0.0     r2_inv*csc(θ)^2
  ]
end


##
# Follows the Formalism of Minkowski Magnification Notes
##

##----------------------------------------------------------------------------------------------------------------------
# Useful functions
##----------------------------------------------------------------------------------------------------------------------

function βboundary(met::Minkowski, α, θo, θs)
  cosθs2 = cos(θs)^2
  √max((cos(θo)^2 - cosθs2) * (α^2) / (cosθs2 - 1), 0.0)
end

"""
  r_potential(met::Minkowski, r, η, λ)

Radial potential of a kerr blackhole

  `met` : Minkowski Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `r`   : Boyer Lindquist radius

"""
function r_potential(met::Minkowski, η, λ, r)
  λ2 = λ^2
  r^2 * (-η - λ2 + r^2) 
end

"""
  θ_potential(met::Minkowski, r, η, λ)

Theta potential of a kerr blackhole

  `met` : Minkowski Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `θ`   : Boyer Lindquist inclination

"""
θ_potential(met::Minkowski, η, λ, θ) = η - λ^2 * cot(θ)^2

"""
  get_roots(b)

Returns double the number `x` plus `1`.
"""
function get_roots(b)
  return @SVector[-b, b]
end

"""
  φ(θs, φ, θo, n)

Screen arival angle of emission.

  `θs`: Emitter Inclination

  `ϕ` : Boyer Lindquist ϕ

  `θo`: Observer Inclination

  'n' : Index of nth image
"""
function φ(θs::Real, ϕ::Real, θo::Real, n::Real)#φ screen coordinate
  p1 = cos(ϕ)sin(θs)
  p2 = cos(θs)sin(θo) + cos(θo)sin(θs)sin(ϕ)
  return (sign(p2) * acos(p1 / √(p1^2 + p2^2))) + n * π
end

function ϕ(φ, θo, θs, isindir::Bool)
  if cos(φ) > 0
    return atan((cos(θo)*tan(φ)-√(cot(θs)^2*sin(θo)^2*(cos(θo)^2-cot(θs)^2*sin(θo)^2+tan(φ)^2)))/(cos(θo)^2-cot(θs)^2*sin(θo)^2)) + (sin(φ) > 0 ? 2*π : 0)
  else
    return atan((cos(θo)*tan(φ)+√(cot(θs)^2*sin(θo)^2*(cos(θo)^2-cot(θs)^2*sin(θo)^2+tan(φ)^2)))/(cos(θo)^2-cot(θs)^2*sin(θo)^2)) + π
  end
end

"""
  ψ(θs, ϕ, θo)

Winding angle of direct emission in the forward raytracing problem.

  `θs`: The half opening angle of the cone that the point lies on.

  `ϕ`: The parameterized angle around the point.

  `θo`: The inclination angle of the point.
"""
function ψ(θs::Real, ϕ::Real, θo::Real, isindir::Bool)#ψ angle swept out 
  cψ = cos(θs) * cos(θo) - sin(θs) * sin(θo) * sin(ϕ)
  sψ = √(cos(ϕ)^2*sin(θs)^2 + (cos(θs)*sin(θo) +cos(θo)*sin(θs)*sin(ϕ))^2)
  return atan(sψ, cψ)
  #return acos(cψ)
end

"""
  ψ(α, β, θs, θo, n)

Winding angle of direct emission in the forward raytracing problem.

  `α`: The half opening angle of the cone that the point lies on.

  `β`: The parameterized angle around the point.

  `θs`: The half opening angle of the cone that the point lies on.

  `θo`: The inclination angle of the point.

  `n`: Image order
"""
function ψ(α, β, θs, θo, isindir::Bool, n)
  signβ = sign(β)
  isincone = abs(cos(θs)) < abs(cos(θo))
  if isincone && (isindir != ((signβ > 0) ⊻ (θo > π / 2)))
    return Inf
  end
  if ((((signβ < 0) ⊻ (θs > π / 2)) ⊻ (n % 2 == 1)) && !isincone )
    return Inf
  end

  φ = atan(β, α)

  νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > π / 2)
  return ψ(θs, ϕ(φ, θo, θs, isindir), θo, isindir) 
  #if angle < 0 && φ > π/2
  #  angle += 2*π
  #end
  return angle
  #Ghat = n*π
  #Go = ψ(θs, atan(sϕ, cϕ), θo) 
  #Gs = 0
  #return real(isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone

end

Δ(r) = r^2
Σ(r, θ) = r^2 * cos(θ)^2
A(r) = r^4 
Ξ(r) = r^4

η(α, β, θo) = (α^2) * cos(θo)^2 + β^2
λ(α, θo) = -α * sin(θo)
##----------------------------------------------------------------------------------------------------------------------
# Radial Stuff
##----------------------------------------------------------------------------------------------------------------------

function rs(met::Minkowski, α, β, θs, o::AssymptoticObserver, isindir::Bool, n::Int64)
  θo = o.inclination
  if cos(θs) > abs(cos(θo))
    βbound = (abs(α) >= eps() ? βboundary(met, α, θo, θs) : 0.0)
    if abs(β) + eps() < βbound
      return 0.0, true, 1
    end
  end
  angle = ψ(α, β, θs, θo, isindir, n)
  if  angle == Inf || angle == -Inf || tan(angle) == Inf 
    return 0.0, true, 1
  end

  return √(α^2 + β^2)*csc(angle), true, 1
end

"""
  rs(b, ψ)

Emission radius of source.

  `b`: Impact parameter

  `ψ`: Winding angle
"""
function rs(b, ψ)
  return b*csc(ψ)
end

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
function p_boyer_lindquist_d(met::Minkowski, r, θ, η, λ, νr::Bool, νθ::Bool)
  @SVector [-1, (νr ? 1 : -1) * √abs(r_potential(met, η, λ, r)) / Δ(r), λ, (νθ ? 1 : -1) * √abs(θ_potential(met, η, λ, θ))]
end

"""
    jac_bl2zamo_du(r, θ, a)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `r` : Radius
    
    `θ` : Inclination 
"""
function jac_bl2zamo_du(met::Minkowski, r, θ)
  # coords = {t, r, ϕ, θ}
  Σt = Σ(r, θ)
  Δt = Δ(r)
  At = A(r)
  return @SMatrix [# Eq 3.1 1972ApJ...178..347B
    √(At / (Σt * Δt)) 0.0         0.0                   0.0
    0.0               √(Δt / Σt)  0.0                   0.0
    0.0               0.0         √(Σt / At)*csc(θ)     0.0
    0.0               0.0         0.0                   -1/√Σt
  ]
end

"""
    jac_zamo2zbl_du(r, θ, a)

Jacobian which converts ZAMO covector on the right to a Boyer-Lindquist (BL) covector

    `r` : Radius
    
    `θ` : Inclination 
"""
function jac_zamo2bl_du(met::Minkowski, r, θ)
  Σt = Σ(r, θ)
  Δt = Δ(r)
  At = A(r)

  return @SMatrix [
    # coords = {t, r, ϕ, θ}
    √((Σt * Δt) / At) 0.0                     0.0                       0.0
    0.0               √(Σt / Δt)              0.0                       0.0
    0.0               0.0                     √(At / Σt)*sin(θ)         0.0
    0.0               0.0                     0.0                       -√Σt
  ]
end


function jac_bl2zamo_ud(met::Minkowski, r, θ)
  Σt = Σ(r, θ)
  Δt = Δ(r)
  At = A(r)

  return @SMatrix [#  Eq 3.2 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
    √((Σt * Δt) / At)             0.0         0.0               0.0
    0.0                           √(Σt / Δt)  0.0               0.0
    0.0                           0.0         √(At / Σt)*sin(θ) 0.0
    0.0                           0.0         0.0               -√Σt
  ]
end

function jac_zamo2bl_ud(met::Minkowski, r, θ)
  Σt = Σ(r, θ)
  Δt = Δ(r)
  At = A(r)


  return @SMatrix [
    # coords = {t, r, ϕ, θ}
    √(At / (Σt * Δt))     0.0         0.0               0.0
    0.0                   √(Δt / Σt)  0.0               0.0
    0.0                   0.0         √(Σt / At)*csc(θ) 0.0
    0.0                   0.0         0.0               -1/√Σt
  ]
end

function jac_zamo2fluid_ud(β, θ, φ)
  γ = 1 / √(1 - β^2)
  sinφ = sin(φ)
  cosφ = cos(φ)
  sinθ = sin(θ)
  cosθ = cos(θ)

  return @SMatrix [
    γ                 -β*γ*cosφ*sinθ                        -β*γ*sinφ*sinθ                    -β*γ*cosθ
    -β*γ*cosφ*sinθ    cosθ^2*cosφ^2+γ*cosφ^2*sinθ^2+sinφ^2  (γ-1)*cosφ*sinθ^2*sinφ            (γ-1)*cosθ*cosφ*sinθ
    -β*γ*sinθ*sinφ    (γ-1)*cosφ*sinθ^2*sinφ                cosφ^2+(cosθ^2+γ*sinθ^2)*sinφ^2   (γ-1)*cosθ*sinθ*sinφ
    -β*γ*cosθ         (γ-1)*cosθ*cosφ*sinθ                  (γ-1)*cosθ*sinθ*sinφ              γ*cosθ^2+sinθ^2
  ]
end

function penrose_walker(met::Minkowski, r, p_u::AbstractVector, f_u::AbstractVector)# Eq 6 arXiv:2001.08750v1
  pt, pr, pϕ, pθ = p_u
  ft, fr, fϕ, fθ = f_u

  A = pt * fr - pr * ft
  B = r^2 * (pϕ * fθ - pθ * fϕ) 
  return A * r, B * r
end

function screen_polarisation(met::Minkowski, κ::Complex, α, β)# Eq 31 10.1103/PhysRevD.104.044060
  #TODO: Check which is real and which is imaginary
  κ1 = real(κ)
  κ2 = imag(κ)

  μ = -α
  fα = (β * κ2 - μ * κ1) / (μ^2 + β^2)
  fβ = (β * κ1 + μ * κ2) / (μ^2 + β^2)


  return fα, fβ
end

evpa(fα, fβ) = atan(-fα, fβ)

function calcPol(met::Minkowski, α, β, ri, θs, θo, cross_spec_index, magfield::SVector{3, Float64}, βfluid::SVector{3, Float64}, νr::Bool, θsign::Bool)
  βv = βfluid[1]
  θz = βfluid[2]
  ϕz = βfluid[3]

  ηtemp = η(α, β, θo)
  λtemp = λ(α, θo)
  p_bl_d = p_boyer_lindquist_d(met, ri, θs, ηtemp, λtemp, νr, θsign)

  p_bl_u = met_uu(met, ri, θs) * p_bl_d
  p_zamo_u = jac_bl2zamo_ud(met, ri, θs) * p_bl_u
  p_fluid_u = jac_zamo2fluid_ud(βv, θz, ϕz) * p_zamo_u
  magfieldx, magfieldy, magfieldz = magfield
  _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
  vec = @SVector[p_fluid_uy * magfieldz - p_fluid_uz * magfieldy, p_fluid_uz * magfieldx - p_fluid_ux * magfieldz, p_fluid_ux * magfieldy - p_fluid_uy * magfieldx]
  norm = √sum(vec .* vec) + eps()
  f_fluid_u = SVector(zero(eltype(vec)), vec[1], vec[2], vec[3])
  f_zamo_u = jac_zamo2fluid_ud(-βv, θz, ϕz) * f_fluid_u
  f_bl_u = jac_zamo2bl_ud(met, ri, θs) * f_zamo_u
  A = @SMatrix [
    0.0   1.0  0.0  0.0
    -1.0  0.0  0.0  0.0
    0.0   0.0  0.0  0.0
    0.0   0.0  0.0  0.0
  ]
  B = @SMatrix [
    0.0   0.0  0.0      0.0 
    0.0   0.0  0.0      0.0
    0.0   0.0  0.0      (ri^2)
    0.0   0.0  -(ri^2)  0.0
  ]

  f_temp_d = ((A - B * im) * ri) * (f_bl_u)
  κ = sum(p_bl_u .* f_temp_d)
  κ = κ / √(conj(κ) * κ)


  eα, eβ = screen_polarisation(met, κ, α, β) .* (norm^((cross_spec_index + 1.0) / 2))

  return eα, eβ, 1 / p_fluid_u[1], abs(p_fluid_u[1] / p_fluid_u[4])
end

