# Metric Definition
##----------------------------------------------------------------------------------------------------------------------
struct Schwarzschild <: AbstractMetric end
horizon(met::Schwarzschild) = 2

"""
    met_uu(r, θ, a)

Inverse Schwarzschild Metric in Boyer Lindquist (BL)/Spherical Polar coordinates.

    `met` : Kerr Metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function met_uu(met::Schwarzschild, r, θ)
  r_inv = 1/r
  f = 1 - 2*r_inv
  r2_inv = r_inv^2

  return @SMatrix [
    1/f   0.0   0.0     0.0
    0.0   f     0.0     0.0
    0.0   0.0   r2_inv  0.0
    0.0   0.0   0.0     r2_inv*csc(θ)^2
  ]
end

"""
    met_dd(r, θ, a)

Schwarzschild Metric in Boyer Lindquist (BL)/Spherical Polar coordinates.

    `met` : Kerr Metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function met_dd(met::Schwarzschild, r, θ)
  r_inv = 1/r
  f = 1 - 2*r_inv
  r2 = r^2

  return @SMatrix [
    -f   0.0   0.0     0.0
    0.0   1/f     0.0     0.0
    0.0   0.0   r2  0.0
    0.0   0.0   0.0     r2*sin(θ)^2
  ]
end



##
# Follows the Formalism of Schwarzschild Magnification Notes
##

##----------------------------------------------------------------------------------------------------------------------
# Useful functions
##----------------------------------------------------------------------------------------------------------------------

function βboundary(met::Schwarzschild, α, o::AssymptoticObserver, θs)
  θo = o.inclination
  cosθs2 = cos(θs)^2
  √max((cos(θo)^2 - cosθs2) * α^2 / (cosθs2 - 1), 0.0)
end

"""
  r_potential(met::Schwarzschild, r, η, λ)

Radial potential of a kerr blackhole

  `met` : Schwarzschild Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `r`   : Boyer Lindquist radius

"""
function r_potential(met::Schwarzschild, η, λ, r)
  λ2 = λ^2
  r * (2 * η + 2 * λ2 + r * (-η - λ2 + r^2)) # Eq 7 PhysRevD.101.044032
end

"""
  θ_potential(met::Schwarzschild, r, η, λ)

Theta potential of a kerr blackhole

  `met` : Schwarzschild Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `θ`   : Boyer Lindquist inclination

"""
θ_potential(met::Schwarzschild, η, λ, θ) = η - λ^2 * cot(θ)^2

"""
  get_roots(b)

Returns double the number `x` plus `1`.
"""
function get_roots(b)
  q = 2 * b^2 + 0im
  p = -b^2 + 0im
  u = (-q / 2 + (q^2 / 4 + p^3 / 27)^(1 / 2))^(1 / 3)
  C = ((-1 + 0im)^(2 / 3), (-1 + 0im)^(4 / 3), 1) .* u
  v = -p ./ ((3 + 0im) .* C)
  return sort([i for i in C .+ v], lt=(x, y) -> real(x) < real(y))
end

"""
  ϕ(θs, φ, θo, n)

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

"""
  p(α::Real, φ::Real, θ::Real, χ::Real)

Parameterized point in space bulk spacetime. The point is assumed to lie on a cone whose apex coincides with the 
center of the blackhole.

  `α`: The half opening angle of the cone that the point lies on.

  `φ`: The parameterized angle around the point.

  `θ`: The inclination angle of the point.

  `χ`: The azimuthal rotation of the point.
"""
function p(α, φ, θ, χ)

  return (
    cos(χ) * cos(φ) * sin(α) + sin(χ) * (cos(α) * sin(θ) - cos(θ) * sin(α) * sin(φ)),
    sin(χ) * cos(φ) * sin(α) - cos(χ) * (cos(α) * sin(θ) - cos(θ) * sin(α) * sin(φ)),
    cos(α) * cos(θ) + sin(α) * sin(θ) * sin(φ)
  )
end

#TODO: Fix this
function ϕ(φ, θo, θs, isindir::Bool)
  
  if cos(φ) > 0
    return atan((cos(θo)*tan(φ)-√(cot(θs)^2*sin(θo)^2*(cos(θo)^2-cot(θs)^2*sin(θo)^2+tan(φ)^2)))/(cos(θo)^2-cot(θs)^2*sin(θo)^2))
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
  cϕ = csc(θs) * cos(φ)
  sϕ = csc(θs) * sec(θo) * (sin(φ) - cos(θs) * sin(θo))

  νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > π / 2)


  #angle = ψ(θs, atan(sϕ, cϕ), θo, isindir) 
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

Δ(r) = r^2 - 2r
Σ(r, θ) = r^2 * cos(θ)^2
A(r) = r^4 
Ξ(r) = r^4


η(α, β, θo) = (α^2) * cos(θo)^2 + β^2
λ(α, θo) = -α * sin(θo)
##----------------------------------------------------------------------------------------------------------------------
# Radial Stuff
##----------------------------------------------------------------------------------------------------------------------

function rs(met::Schwarzschild, α, β, θs, o::AssymptoticObserver, isindir::Bool, n::Int64)
  θo = o.inclination
  if cos(θs) > abs(cos(θo))
    βbound = (abs(α) >= eps() ? βboundary(met, α, o, θs) : 0.0)
    if abs(β) + eps() < βbound
      return 0.0, true, 4
    end
  end
  angle = ψ(α, β, θs, θo, isindir, n)
  if  angle == Inf || angle == -Inf || tan(angle) == Inf 
    return 0.0, true, 4
  end

  #return √(α^2 + β^2)*csc(angle)
  winding_angle = ψ(α, β, θs, θo, isindir, n) 
  #TODO: Figure out what this should be. The signature is rs, νr, number of real roots
  return rs(√(α^2 + β^2), winding_angle), true, 4 
end

"""
  rs(b, ψ)

Emission radius of source.

  `b`: Impact parameter

  `ψ`: Winding angle
"""
function rs(b, ψ)
  rs = b^2 > 27 ? _rsout(b, ψ) : _rsin(b, ψ)
  return rs < 2 ? 0.0 : rs
end

"""
  _rsout(b, ψ)

Emission radius for emission that lies outside the photon ring
"""
function _rsout(b, ψ)
  r1, r3, r4 = get_roots(b)
  r32 = r3
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4
  k = abs((r32 * r41) / (r31 * r42))
  fo = FastElliptic.F(asin(real(√(r31 / r41))), min(k, 1))
  mino_time = ψ / b
  arg = fo - real(1 / 2 * √(r31 * r42) * mino_time)
  arg ≤ -fo + eps() && return 0.0
  sn = r41 * FastElliptic.sn(arg, k)^2
  return real((r31 * r4 - r3 * sn) / (r31 - sn))
end

"""
  _rsin(b, ψ)

Emission radius for emission that lies inside the photon ring

  `b` : Impact parameter

  `ψ` : Winding angle
"""
function _rsin(b, ψ)
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r21 = r2 - r1
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4
  A = √(r32 * r42)
  B = √(r31 * r41)
  k = min(real(((A + B)^2 - r21^2) / (4 * A * B)), 1)
  fo = FastElliptic.F(real(acos((A - B) / (A + B))), k)
  mino_time = ψ / b
  arg = real(fo - √(A * B) * mino_time)
  arg ≤ - eps() && return 0.0
  cn = FastElliptic.cn(arg, k)
  num = -A * r1 + B * r2 + (A * r1 + B * r2) * cn
  den = -A + B + (A + B) * cn
  return real(num / den)
end

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
function p_boyer_lindquist_d(met::Schwarzschild, r, θ, η, λ, νr::Bool, νθ::Bool)
  @SVector [-1, (νr ? 1 : -1) * √abs(r_potential(met, η, λ, r)) / Δ(r), λ, (νθ ? 1 : -1) * √abs(θ_potential(met, η, λ, θ))]
end


"""
    jac_bl2zamo_du(met::Schwarzschild, r, θ)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `met` : Kerr Metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function jac_bl2zamo_du(met::Schwarzschild, r, θ)
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
    jac_zamo2zbl_du(met::Schwarzschild, r, θ)

Jacobian which converts ZAMO covector on the right to a Boyer-Lindquist (BL) covector

    `met` : Kerr Metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function jac_zamo2bl_du(met::Schwarzschild, r, θ)
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


function jac_bl2zamo_ud(met::Schwarzschild, r, θ)
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

function jac_zamo2bl_ud(met::Schwarzschild, r, θ)
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

function penrose_walker(met::Schwarzschild, r, p_u::AbstractVector, f_u::AbstractVector)# Eq 6 arXiv:2001.08750v1
  pt, pr, pϕ, pθ = p_u
  ft, fr, fϕ, fθ = f_u

  A = pt * fr - pr * ft
  B = r^2 * (pϕ * fθ - pθ * fϕ) 
  return A * r, B * r
end

function screen_polarisation(met::Schwarzschild, κ::Complex, α, β)# Eq 31 10.1103/PhysRevD.104.044060
  #TODO: Check which is real and which is imaginary
  κ1 = real(κ)
  κ2 = imag(κ)

  μ = -α
  fα = (β * κ2 - μ * κ1) / (μ^2 + β^2)
  fβ = (β * κ1 + μ * κ2) / (μ^2 + β^2)


  return fα, fβ
end

evpa(fα, fβ) = atan(-fα, fβ)

function calcPol(met::Schwarzschild, α, β, ri, θs, θo, cross_spec_index, magfield::SVector{3, Float64}, βfluid::SVector{3, Float64}, νr::Bool, θsign::Bool)
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

##----------------------------------------------------------------------------------------------------------------------
# Magnification Stuff
##----------------------------------------------------------------------------------------------------------------------

"""
  g(rs::real)

Redshift factor
  `rs` : The emission radius of the photon
"""
g(rs::Real) = √(1 - 2 / rs)

dψ_drdb(b::Real, r::Real)::Real = real(r^4 / ((r^4 - b^2 * r * (r - 2 + 0im))^(3 / 2)))

function ψt(b::Real)
  if b^2 < 27
    return Inf
  else
    r1, r3, r4 = get_roots(b)
    r2 = 0
    r32 = r3 - r2
    r41 = r4 - r1
    r31 = r3 - r1
    r42 = r4 - r2

    k = min(real((r32 * r41 / (r42 * r31))), 1)
    g = 2 / √(r42 * r31)
    return real(b * g * Elliptic.F(real(asin(√(r31 / r41))), k))
  end
end

function dψ_du(b::Real, u::Real)
  """
      u^2 = r-r4
  """
  r1, r3, r4 = get_roots(b)
  r2 = 0

  return 2 * b / √((u^2 + r4 - r1) * (u^2 + r4 - r2) * (u^2 + r4 - r3))
end

function droots_db(b::Real)
  p = -b^2
  q = 2 * b^2
  disc = (q^2 / (4 + 0im) + p^3 / (27 + 0im))^(1 / 2)
  C = (-q / (2 + 0im) + disc)^(1 / 3)

  temp_desc = √(81 + 0im - 3 * b^2)
  dC_db = -(3)^(1 / 3) * (2 * b + b * (b^2 - 18) / temp_desc) / (b^2 * (temp_desc - 9))^(2 / 3)
  dp_db = -2 * b
  dp_3C_db = (dp_db / C - p * dC_db / C^2) / 3

  dr1_db = dC_db * (-1 + √(3) * 1im) / 2 - dp_3C_db / ((-1 + √(3) * 1im) / 2)
  dr3_db = dC_db * (-1 - √(3) * 1im) / 2 - dp_3C_db / ((-1 - √(3) * 1im) / 2)
  dr4_db = dC_db - dp_3C_db
  return dr1_db, dr3_db, dr4_db
end

"""
  ψbackwards_out_front(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem on the front side of a ring

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_out_front(b::Real, rs::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  k = min(real((r32 * r41) / (r31 * r42)), 1)
  g = 2 / √(r42 * r31)
  fo = Elliptic.F(real(asin(√(r31 / r41))), k)
  return real(b * real(g * (fo - Elliptic.F(real(asin(√(r31 * (rs - r4) / (r41 * (rs - r3))))), k))))
end

"""
  ψbackwards_out_bank(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem on the back side of a ring

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_out_back(b::Real, rs::Real)::Real
  return real(2ψt(b) - ψbackwards_out_front(b, rs))
end

"""
  ψbackwards_in(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem for rays that fall into the blackhole

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_in(b::Real, rs::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2
  r21 = r2 - r1

  A = √(r32 * r42)
  B = √(r31 * r41)
  k = min(real(((A + B)^2 - r21^2) / (4 * A * B)), 1)
  fo = Elliptic.F(real(acos((A - B) / (A + B))), k)
  x = real((A * (rs - r1) - B * (rs - r2)) / (A * (rs - r1) + B * (rs - r2)))
  return real(b * (fo - Elliptic.F(acos(x), k)) / √(A * B))
end


dsnInv_dx(x::Real, k::Real) = 1 / √((1 - x^2) * (1 - k * x^2))
dsnInv_dk(x::Real, k::Real) = -(2Elliptic.E(asin(x), k) + 2(k - 1)Elliptic.F(asin(x), k) - k * sin(2asin(x)) / √(1 - k * x^2)) / (4(k - 1) * k)
dcnInv_dx(x::Real, k::Real) = -1 / √((1 - x^2) * (1 - k * (1 - x^2)))
dcnInv_dk(x::Real, k::Real) = -(2Elliptic.E(acos(x), k) + 2(k - 1)Elliptic.F(acos(x), k) - k * sin(2acos(x)) / √(1 + k * (x^2 - 1))) / (4(k - 1) * k)

function dgfo_out_db(b::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)

  k = min(real(r32 * r41 / (r42 * r31)), 1)
  dk = (dr3 / r32 + (dr4 - dr1) / r41 - dr4 / r42 - (dr3 - dr1) / r31) * k
  g = 2 / √(r42 * r31)
  dg = -(dr4 / r42 + (dr3 - dr1) / r31) * g / 2
  x = real(√(r31 / r41))
  dx = ((dr3 - dr1) / r31 - (dr4 - dr1) / r41) * x / 2

  ans = dg * Elliptic.F(real(asin(x)), k)
  ans += g * dsnInv_dx(x, k) * dx
  ans += g * dsnInv_dk(x, k) * dk
  return real(ans)
end

function dΔτout_db(b::Real, rs::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)


  k = min(real(r32 * r41 / (r42 * r31)), 1)
  dk_db = (dr3 / r32 + (dr4 - dr1) / r41 - dr4 / r42 - (dr3 - dr1) / r31) * k
  g = 2 / √(r42 * r31)
  dg = -(dr4 / r42 + (dr3 - dr1) / r31) * (g / 2)
  x = real(√(r31 * (rs - r4) / (r41 * (rs - r3))))
  dx_db = ((dr3 - dr1) / r31 - dr4 / (rs - r4) - (dr4 - dr1) / r41 + dr3 / (rs - r3)) * x / 2

  ans = dg * Elliptic.F(asin(x), k)
  ans += g * dsnInv_dx(x, k) * dx_db
  ans += g * dsnInv_dk(x, k) * dk_db

  return real(ans)
end

function dgfo_in_db(b::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r21 = r2 - r1
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)
  dr2 = 0

  A = √(r32 * r42)
  dA = (dr3 / r32 + dr4 / r42) * (A / 2)
  B = √(r31 * r41)
  dB = ((dr3 - dr1) / r31 + (dr4 - dr1) / r41) * (B / 2)
  g = 1 / √(A * B)
  dg = (-dA / A - dB / B) * (g / 2)
  x = real((A - B) / (A + B))
  dx = ((dA - dB) - (A - B) * (dA + dB) / (A + B)) / (A + B)
  k = min(real(((A + B)^2 - r21^2) / (4A * B)), 1)
  dk = (2((dA + dB) * (A + B) - (dr2 - dr1) * (r21)) + ((A + B)^2 - r21^2) * (-dA / A - dB / B)) / (4A * B)

  ans = dg * Elliptic.F(acos(x), k)
  ans += g * dcnInv_dx(x, k) * dx
  ans += g * dcnInv_dk(x, k) * dk
  return real(ans)
end


function dΔτin_db(b::Real, rs::Real)
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2
  r21 = r2 - r1

  dr1, dr3, dr4 = droots_db(b)
  dr2 = 0

  A = √(r32 * r42)
  dA = (dr3 / r32 + dr4 / r42) * (A / 2)
  B = √(r31 * r41)
  dB = ((dr3 - dr1) / r31 + (dr4 - dr1) / r41) * (B / 2)
  k = min(real(((A + B)^2 - r21^2) / (4A * B)), 1)
  dk = (2((dA + dB) * (A + B) - (dr2 - dr1) * (r21)) + ((A + B)^2 - r21^2) * (-dA / A - dB / B)) / (4A * B)
  g = 1 / √(A * B)
  dg = (-dA / A - dB / B) * (g / 2)
  x = real((A * (rs - r1) - B * (rs - r2)) / (A * (rs - r1) + B * (rs - r2)))
  dx = ((dA * (rs - r1) - A * dr1 - dB * (rs - r2)) - (A * (rs - r1) - B * (rs - r2)) * (dA * (rs - r1) - A * dr1 + dB * (rs - r2)) / (A * (rs - r1) + B * (rs - r2))) / (A * (rs - r1) + B * (rs - r2))

  ans = dg * Elliptic.F(acos(x), k)
  ans += g * dcnInv_dx(x, k) * dx
  ans += g * dcnInv_dk(x, k) * dk
  return real(ans)
end

dψout_front_db(b::Real, rs::Real) = b * (dgfo_out_db(b) - dΔτout_db(b, rs)) + ψbackwards_out_front(b, rs) / b
dψout_back_db(b::Real, rs::Real) = b * (dgfo_out_db(b) + dΔτout_db(b, rs)) + ψbackwards_out_back(b, rs) / b
dψin_db(b::Real, rs::Real) = b * (dgfo_in_db(b) - dΔτin_db(b, rs)) + ψbackwards_in(b, rs) / b


function dψ_dudb(b::Real, u::Real)
  """
      u^2 = r-r4
  """
  r1, r3, r4 = get_roots(b)
  r2 = 0
  dr1_db, dr3_db, dr4_db = droots_db(b)
  dr2_db = 0
  n1 = -1 / 2 * (dr4_db - dr1_db) / (u^2 + r4 - r1)
  n2 = -1 / 2 * (dr4_db - dr2_db) / (u^2 + r4 - r2)
  n3 = -1 / 2 * (dr4_db - dr3_db) / (u^2 + r4 - r3)

  d = 2 / √((u^2 + r4 - r1) * (u^2 + r4 - r2) * (u^2 + r4 - r3))

  return d * (b * (n1 + n2 + n3) + 1)
end

"""
  magnification(b::Real, ψ::Real)::Real

Analytic value for the magnifications of a point lensed by a blackhole

  `b` : Real

  `ψ` : Real
"""
function magnification(b::Real, ψ::Real)::Real
  ri = rs(b, ψ)
  dψ_db = 0
  if b^2 >= 27
    if ψ > ψt(b)
      dψ_db = dψout_back_db(b, ri)
    else
      dψ_db = dψout_front_db(b, ri)
    end
  else
    dψ_db = dψin_db(b, ri)
  end
  return abs(b / √radial_potential(b, ri) / (sin(ψ) * dψ_db))
end

"""
  magnification_numeric(b::Real, ψ::Real)::Real

Numeric value for the magnifications of a point lensed by a blackhole

  `b` : Real

  `ψ` : Real
"""
function magnification_numeric(b::Real, ψ::Real)::Real
  ri = rs(b, ψ)
  dψ_db = quadgk(r -> dψ_drdb(b, r), ri, 1000)[1]
  if b^2 > 27 && ψ > ψt(b)
    p = -b^2
    q = 2b^2
    disc = √(q^2 / 4 + p^3 / 27 + 0im)
    C = (-q / 2 + disc)^(1 / 3)

    temp_desc = √(81 - 3b^2 + 0im)
    dC_db = -3^(1 / 3) * (2b + b * (b^2 - 18) / temp_desc) / (b^2 * (temp_desc - 9))^(2 / 3)
    dp_db = -2 * b
    dp_3C_db = (dp_db / C - p * dC_db / C^2) / 3

    rT = C - p / (3 * C)
    us = real(√(ri - rT))
    dψ_db += 2 * quadgk(u -> dψ_dudb(b, real(u)), 0, us)[1]
    dψ_db -= dψ_du(b, us) / us * (dC_db - dp_3C_db)
  end
  return abs(b / √radial_potential(b, ri) / (sin(ψ) * dψ_db))
end

ψequatorial(φ::Real, θ::Real) = acos(-sin(θ) * sin(φ))
ϕequatorial(φ::Real, θ::Real) = sign(sin(φ)) * acos(cos(φ) / (√((cos(θ) * sin(φ))^2 + cos(φ)^2)))
