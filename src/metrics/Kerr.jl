# metric Definition
##----------------------------------------------------------------------------------------------------------------------
function elliptic_pi(n::Float64, φ::Float64, m::Float64)
  return ArbNumerics.elliptic_pi(ArbFloat(n), ArbFloat(φ), ArbFloat(m))
end
function elliptic_pi(n::Float64, m::Float64)
  return ArbNumerics.elliptic_pi(ArbFloat(n), ArbFloat(m))
end

struct Kerr <: AbstractMetric
    spin::Float64
end

"""
    horizon(metric::Kerr)

Horizon for the Kerr metric.

    `metric` : Kerr metric

"""
horizon(metric::Kerr) = 1 + √(1 - metric.spin^2)

Δ(metric::Kerr, r, a) = r^2 - 2r + a^2
Σ(metric::Kerr, r, θ, a) = r^2 + a^2 * cos(θ)^2
A(metric::Kerr, r, θ, a) = (r^2 + a^2)^2 - a^2 * Δ(metric, r, a) * sin(θ)^2
Ξ(metric::Kerr, r, θ, a) = (r^2 + a^2)^2 - Δ(metric, r, a) * a^2 * sin(θ)^2
ω(metric::Kerr, r, θ, a) = 2 * a * r / Ξ(metric, r, θ, a)

η(metric::Kerr, α, β, θo, a) = (α^2 - a^2) * cos(θo)^2 + β^2
λ(metric::Kerr, α, θo) = -α * sin(θo)

rtildep(a) = 2 * (1 + Cos(2 / 3 * acos(a)))
rtilden(a) = 2 * (1 + Cos(2 / 3 * acos(-a)))

"""
    metric_uu(metric::Kerr, r, θ)

Inverse Kerr metric in Boyer Lindquist (BL) coordinates.

    `metric` : Kerr metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function metric_uu(metric::Kerr, r, θ)
  a = metric.spin
  Ξt = Ξ(metric, r, θ, a)
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  ωt = ω(metric, r, θ, a)

  return @SMatrix [ #Eq 1 2105.09440
    -Ξt/(Σt*Δt)     0.0     -Ξt*ωt/(Σt*Δt)                  0.0
    0.0             Δt/Σt   0.0                             0.0
    -Ξt*ωt/(Σt*Δt)  0.0     Σt*csc(θ)^2/Ξt-Ξt*ωt^2/(Σt*Δt)  0.0
    0.0             0.0     0.0                             1/Σt
  ]
end

"""
    metric_dd(metric::Kerr, r, θ)

Inverse Kerr metric in Boyer Lindquist (BL) coordinates.

    `metric` : Kerr metric

    `r` : Radius
    
    `θ` : Inclination 
"""
function metric_dd(metric::Kerr, r, θ)
  a = metric.spin
  Ξt = Ξ(metric, r, θ, a)
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  ωt = ω(metric, r, θ, a)
  sin2 = sin(θ)^2

  return @SMatrix [ #Eq 1 2105.09440
    -(Δt*Σt/Ξt)+(Ξt*ωt^2*sin2/Σt)     0.0       -Ξt*ωt*sin2/Σt     0.0
    0.0                               Σt/Δt     0.0                0.0
    -Ξt*ωt*sin2/Σt                    0.0       Ξt*sin2/Σt         0.0
    0.0                               0.0       0.0                Σt
  ]
end

# Follows the Formalism of Gralla & Lupsasca (https://arxiv.org/pdf/1910.12881.pdf)

# Useful functions
##----------------------------------------------------------------------------------------------------------------------

αboundary(metric::Kerr, θs) = metric.spin * sin(θs)

"""
    βboundary(metric::Kerr, α, o::AssymptoticObserver, θs)

Defines a boundary on the Assyptotic observers screen where emission that originates from θs must fall within.

    `metric`: Kerr metric

    `α`   : Bardeen horizontal coordinate

    `o`  : Assymptotic Observer

    `r`   : Boyer Lindquist radius

"""
function βboundary(metric::Kerr, α, o::AssymptoticObserver, θs) 
  θo = o.inclination
  a = metric.spin
  cosθs2 = cos(θs)^2
  √max((cos(θo)^2 - cosθs2) * (α^2 - a^2*(1 - cosθs2)) / (cosθs2 - 1), 0.0) #eq 15 DOI 10.3847/1538-4357/acafe3 
end

"""
    r_potential(metric::Kerr, r, η, λ)

Radial potential of spacetime

    `metric`: Kerr metric

    `η`  : Reduced Carter constant

    `λ`  : Reduced azimuthal agular momentum

    `r`  : Boyer Lindquist radius

"""
function r_potential(metric::Kerr, η, λ, r)
  a = metric.spin
  λ2 = λ^2
  a * (a * (r * (r + 2) - η) - 4 * λ * r) + r * (2 * η + 2 * λ2 + r * (-η - λ2 + r^2)) # Eq 7 PhysRevD.101.044032
end

"""
    θ_potential(metric::Kerr, r, η, λ)

Theta potential of a kerr blackhole

    `metric`: Kerr metric

    `η`  : Reduced Carter constant

    `λ`  : Reduced azimuthal agular momentum

    `θ`  : Boyer Lindquist inclination

"""
function θ_potential(metric::Kerr, η, λ, θ)
  a = metric.spin
  η + a^2 * cos(θ)^2 - λ^2 * cot(θ)^2
end

"""
  get_radial_roots(metric::Kerr, η, λ)

Returns roots of r⁴ + (a²-η-λ²)r² + 2(η+(a-λ)²)r - a²η

  `metric`: Kerr metric

  `η`  : Reduced Carter constant

  `λ`  : Reduced azimuthal agular momentum

"""
function get_radial_roots(metric::Kerr, η, λ)
  a = metric.spin

  a2 = a * a
  A = a2 - η - λ * λ
  B = 2(η + (λ - a)^2)
  C = -a2 * η

  P = -A * A / 12 - C
  Q = -A / 3 * (A * A / 36 + 0im - C) - B * B / 8

  Δ3 = -4 * P^3 - 27 * Q^2
  ωp = (-Q / 2 + sqrt(-Δ3 / 108) + 0im)^(1 / 3)

  #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
  C = (-0.4999999999999999 + 0.8660254037844387im, -0.5000000000000002 - 0.8660254037844385im, 1) .* ωp

  v = -P ./ ((3 + 0im) .* C)

  ξ0 = argmax(real, (C .+ v)) - A / 3.0

  r1 = (-sqrt(2ξ0) - sqrt(-(2A + 2ξ0 - (√2B) / (√ξ0)))) / 2
  r2 = (-sqrt(2ξ0) + sqrt(-(2A + 2ξ0 - (√2B) / (√ξ0)))) / 2
  r3 = (sqrt(2ξ0) - sqrt(-(2A + 2ξ0 + (√2B) / (√ξ0)))) / 2
  r4 = (sqrt(2ξ0) + sqrt(-(2A + 2ξ0 + (√2B) / (√ξ0)))) / 2

  return r1, r2, r3, r4
end


"""
    λcrit(r::Complex, a)

Returns λ values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
λcrit(r, a) = a + r / a * (r - 2Δ(metric, r, a) / (r - 1))
"""
    ηcrit(r::Complex, a)

Returns η values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
ηcrit(r, a) = (r^3 / a^2) * (4 * Δ(metric, r, a) / (r - 1)^2 - r)

##----------------------------------------------------------------------------------------------------------------------
# Radial Stuff
##----------------------------------------------------------------------------------------------------------------------

get_root_diffs(r1, r2, r3, r4) = r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2

function rs(metric::Kerr, α, β, θs, o::AssymptoticObserver, isindir, n)
  a = metric.spin
  θo = o.inclination

  if cos(θs) > abs(cos(θo))
    αmin = αboundary(metric, θs)
    βbound = (abs(α) >= αmin + eps() ? βboundary(metric, α, o, θs) : 0.0)
    if abs(β) + eps() < βbound
      #return 0.0, true, 4
      return 0.0
    end
  end
  ηtemp = η(metric, α, β, θo, a)
  λtemp = λ(metric, α, θo)
  τ, _, _ = _Gθ(metric, sign(β), θs, θo, isindir, n, ηtemp, λtemp)
  if τ != Inf
    return _rs(metric, ηtemp, λtemp, τ)[1]
  else
    #return 0.0
    return (0.0, true, 4)
  end
end

function rs_mask(metric::Kerr, n_init, α, β, θs, o::AssymptoticObserver, isindir)
  θo = o.inclination
  a = metric.spin
  if abs(cos(θs)) > abs(cos(θo))
    αmin = αboundary(metric, θs)
    βbound = (abs(α) >= αmin ? βboundary(metric, α, o, θs) : 0.0)
    abs(β) < βbound && return (0.0, true, 4), true
  end
  ηtemp = η(metric, α, β, θo, a)
  λtemp = λ(metric, α, θo)
  #τ = Gθ(α, β, a, θs, θo, isindir, n)[1]
  τ0, τhat, _ = _Gθ(metric, sign(β), θs, θo, isindir, n_init, ηtemp, λtemp)
  if τ0 != Inf
    τ = τ0 + τhat
    return _rs_mask(metric, ηtemp, λtemp, τ0, τ)
  else
    return (0.0, true, 4), true
  end
end

function _rs_mask(metric::Kerr, η, λ, τ0, τ)
  a = metric.spin
  ans = 0.0
  ansmask = 0.0
  νr = true
  νrmask = true

  roots = get_radial_roots(metric, η, λ)
  rh = 1 + √(1 - a * a)
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)

  if numreals == 4 #case 1 & 2
    ans, νr = _rs_case1(metric, real.(roots), rh, τ0)
    ansmask, νrmask = _rs_case1(metric, real.(roots), rh, τ)

  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
      roots = (roots[1], roots[4], roots[2], roots[3])
    end

    ans, νr = _rs_case3(metric, roots, τ0)
    ansmask, νrmask = _rs_case3(metric, roots, τ)
  else #case 4
    ans, νr = _rs_case4(metric, roots, rh, τ0)
    ansmask, νrmask = _rs_case4(metric, roots, rh, τ)
  end
  return (ans, νr, numreals), (ansmask <= rh)
end


"""
  _rs(metric::Kerrk, η, λ, τ)

Emission radius for emission that lies outside the photon ring and whose ray intersects the equatorial plane

  `metric`   : Kerr metric

  `η`     : Normalized Carter constant

  `λ`     : Normalized angular momentum

  `τ`     : Mino Time
"""
function _rs(metric::Kerr, η, λ, τ)
  a = metric.spin
  ans = 0.0
  νr = true

  roots = get_radial_roots(metric, η, λ)
  rh = 1 + √(1 - a^2)
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)

  if numreals == 4 #case 1 & 2
    ans, νr = _rs_case1(metric, real.(roots), rh, τ)
  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
      roots = (roots[1], roots[4], roots[2], roots[3])
    end

    ans, νr = _rs_case3(metric, roots, τ)
  else #case 4
    ans, νr = _rs_case4(metric, roots, rh, τ)
  end
  return ans, νr, numreals
end

function _rs_case1(metric::Kerr, roots, rh, τ)
  _, _, r3, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs

  r4 >= rh && τ > 2I2r_turn(metric, root_diffs) && return 0.0, true # invalid case1
  r4 < rh && τ > I2r(metric, roots, root_diffs, rh, true) && return 0.0, true# invalid case2

  k = (r32 * r41) / (r31 * r42)
  fo = FastElliptic.F(asin(√(r31 / r41)), k)
  X2 = fo - √(r31 * r42) * τ / 2
  sn = r41 * FastElliptic.sn(X2, k)^2
  return (r31 * r4 - r3 * sn) / (r31 - sn), X2 > 0
end

function _rs_case3(metric::Kerr, roots, τ)
  r1, r2, _, _ = roots
  root_diffs = get_root_diffs(roots...)
  r21, r31, r32, r41, r42 = root_diffs

  τ > I3r_full(metric, root_diffs) && return 0.0, true

  A = √abs(r32 * r42)
  B = √abs(r31 * r41)
  k = real(((A + B)^2 - r21^2) / (4 * A * B))

  fo = FastElliptic.F(acos((A - B) / (A + B)), k)
  X3 = real(fo - √(A * B) * τ)
  cn = FastElliptic.cn(X3, k)
  num = -A * r1 + B * r2 + (A * r1 + B * r2) * cn
  den = -A + B + (A + B) * cn

  return real(num / den), X3 > 0
end
function _rs_case4(metric::Kerr, roots, rh, τ)
  r1, _, _, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs

  τ > I4r(metric, roots, root_diffs, rh) && return 0.0, true
  a2 = abs(imag(r1))
  b1 = real(r4)
  C = √real(r31 * r42)
  D = √real(r32 * r41)
  k4 = 4 * C * D / (C + D)^2

  go = √(4a2^2 - (C - D)^2) / ((C + D)^2 - 4a2^2)
  fo = 2 / (C + D) * FastElliptic.F(π / 2 + atan(go), k4)
  X4 = (C + D) / 2 * (fo - τ)
  num = go - FastElliptic.sc(X4, k4)
  den = 1 + go * FastElliptic.sc(X4, k4)

  return -(a2 * num / den + b1), X4 > 0
end

#TODO: define radial minotime

function Ir(metric::Kerr, νr::Bool, θo, rs, α, β)
  a = metric.spin
  return Ir(metric::Kerr, νr::Bool, rs, η(metric, α, β, θo, a), λ(metric, α, θo))
end

function Ir(metric::Kerr, νr::Bool, rs, η, λ)
  ans = 0.0

  roots = get_radial_roots(metric, η, λ)
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)
  root_diffs = get_root_diffs(roots...)

  if numreals == 4 #case 1 & 2
    ans = I2r(metric, real.(roots), real.(root_diffs), rs, νr)
  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
      roots = (roots[1], roots[4], roots[2], roots[3])
    end
    ans = I3r(metric, real.(roots), root_diffs, rs)
  else #case 4
    ans = I4r(metric, real.(roots), root_diffs, rs)
  end
  return ans
end

function I2r_turn(metric::Kerr, root_diffs::NTuple{5})
  _, r31, r32, r41, r42 = root_diffs
  k = r32 * r41 / (r31 * r42)
  return 2 / √real(r31 * r42) * FastElliptic.F(asin(√(r31 / r41)), k)
end

function I2r(metric::Kerr, roots::NTuple{4}, root_diffs::NTuple{5}, rs, isindir)
  _, _, r3, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  k = r32 * r41 / (r31 * r42)
  x2_s = √((rs - r4) / (rs - r3) * r31 / r41)
  !(-1 < x2_s < 1) && return 0.0

  Ir_s = 2 / √real(r31 * r42) * FastElliptic.F(asin(x2_s), k)
  Ir_turn = I2r_turn(metric, root_diffs)

  if isindir
    return Ir_turn - Ir_s
  else
    return Ir_turn + Ir_s
  end
end

function I3r_full(metric::Kerr, root_diffs)
  r21, r31, r32, r41, r42 = map(abs, root_diffs)
  A2 = r32 * r42
  B2 = r31 * r41
  if A2 < 0.0 || B2 < 0
    return Inf
  end

  A, B = √A2, √B2
  k3 = ((A + B)^2 - r21^2) / (4A * B)

  temprat = B / A
  x3_turn = real(((A - B) / (A + B)))
  return 1 / √real(A * B) * FastElliptic.F(acos(x3_turn), k3)
end

function I3r(metric::Kerr, roots, root_diffs, rs)
  r1, r2, _, _ = roots
  r21, r31, r32, r41, r42 = root_diffs

  A2 = real(r32 * r42)
  B2 = real(r31 * r41)
  if A2 < 0.0 || B2 < 0
    return 0
  end
  A, B = √A2, √B2

  k3 = real(((A + B)^2 - r21^2) / (4*A * B))
  temprat = B * (rs - r2) / (A * (rs - r1))
  x3_s = real(((1 + 0im - temprat) / (1 + temprat)))
  Ir_s = 1 / √real(A * B) * FastElliptic.F((acos(x3_s)), k3)
  Ir_full = I3r_full(metric, root_diffs)

  return Ir_full - Ir_s
end

function I4r_full(metric::Kerr, roots, root_diffs)
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  arg1 = real(r31 * r42)
  arg2 = real(r32 * r41)
  if arg1 < 0 || arg2 < 0
    return 0
  end
  C = √arg1
  D = √arg2
  k4 = 4C * D / (C + D)^2
  a2 = abs(imag(r1))

  k4 = 4 * C * D / (C + D)^2

  go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - 4a2^2), 0.0)
  return 1 / (C + D) * FastElliptic.F(π / 2 + atan(go), k4)


end

function I4r(metric::Kerr, roots, root_diffs, rs)
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  if real(r32 * r41) < 0 || real(r31 * r42) < 0
    return 0
  end
  C = √real(r31 * r42)
  D = √real(r32 * r41)
  k4 = 4C * D / (C + D)^2
  a2 = abs(imag(r1))
  b1 = real(r4)

  k4 = 4 * C * D / (C + D)^2

  go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - 4a2^2), 0.0)
  x4_s = (rs + b1) / a2
  Ir_s = 2 / (C + D) * FastElliptic.F(atan(x4_s) + atan(go), k4)
  Ir_full = I4r_full(metric, roots, root_diffs)

  return Ir_full - Ir_s
end

function I2ϕ(metric::Kerr, roots::NTuple{4}, root_diffs::NTuple{5}, rs, isindir, λ)
  _, _, r3, r4 = roots
  r21, r31, r32, r41, r42 = root_diffs
  r43 = r4-r3
  a = metric.spin
  rp = 1 + √(1-a^2)
  rm = 1 - √(1-a^2)
  rp3 = rp - r3
  rp4 = rp - r4
  rm3 = rm - r3
  rm4 = rm - r4

  k = r32 * r41 / (r31 * r42)
  x2_s = √((rs - r4) / (rs - r3) * r31 / r41)
  x2_o = √(r31 / r41)
  !(-1 < x2_s < 1) && return 0.0

  Fo = 2 / √real(r31 * r42) * FastElliptic.F(asin(x2_o), k)
  Fs = 2 / √real(r31 * r42) * FastElliptic.F(asin(x2_s), k)


  coef_p = 2/√(r31*r42)*r43/(rp3*rp4)
  coef_m = 2/√(r31*r42)*r43/(rm3*rm4)
  Πp_s = coef_p*Elliptic.Pi(rp3*r41/(rp4*r31), asin(x2_s), k)
  Πp_turn = coef_p*Elliptic.Pi(rp3*r41/(rp4*r31), asin(x2_o), k)
  Πm_s = coef_m*Elliptic.Pi(rm3*r41/(rm4*r31), asin(x2_s), k)
  Πm_turn = coef_m*Elliptic.Pi(rm3*r41/(rm4*r31), asin(x2_o), k)

  Ip_s = -Πp_s - Fs/rp3
  Ip_turn = -Πp_turn - Fo/rp3
  Im_s = -Πm_s - Fs/rm3
  Im_turn = -Πm_turn - Fo/rm3

  Ip = Ip_turn #- Ip_s
  Im = Im_turn #- Im_s
  #if !isindir
    Ip -= Ip_s
    Im -= Im_s
  #else
  #  Ip += Ip_s
  #  Im += Im_s
  #end

  return 2a/(rp-rm)*((rp-a*λ/2)*Ip - (rm-a*λ/2)*Im)
end

function p1(α, j)
  √((α^2-1)/(j+(1-j)*α^2))
end

function f1(α, sinφ, j)
  p1temp = p1(α, j)
  tempsinφ = √(1-j*sinφ^2)
  return p1temp/2*log(abs((p1temp*tempsinφ + sinφ)/(p1temp*tempsinφ - sinφ)))
end

function R1(α, φ, j)
  return 1/(1-α^2)*(elliptic_pi(α^2/(α^2-1), φ, j) - α*f1(α, sin(φ), j))
end

function I3ϕ(metric::Kerr, roots::NTuple{4}, root_diffs::NTuple{5}, rs, isindir, λ)
  r1, r2, _, _ = roots
  r21, r31, r32, r41, r42 = root_diffs
  a = metric.spin
  rp = 1 + √(1-a^2)
  rm = 1 - √(1-a^2)
  rp1 = real(rp - r1)
  rp2 = real(rp - r2)
  rm1 = real(rm - r1)
  rm2 = real(rm - r2)

  A2 = real(r32 * r42)
  B2 = real(r31 * r41)
  if A2 < 0.0 || B2 < 0
    return 0
  end
  A, B = √A2, √B2

  k3 = real(((A + B)^2 - r21^2) / (4*A * B))
  temprat = B * (rs - r2) / (A * (rs - r1))
  x3_s = real(((1 - temprat) / (1 + temprat)))
  x3_o = (A-B)/(A+B)
  φ_s = acos(x3_s)
  φ_o = acos(x3_o)

  αp = (B*rp2+A*rp1)/(B*rp2-A*rp1)
  αm = (B*rm2+A*rm1)/(B*rm2-A*rm1)

  Fr_s = 1 / √real(A * B) * FastElliptic.F(φ_s, k3)
  Fr_turn = 1 / √real(A * B) * FastElliptic.F(φ_o, k3)

  R1p_s = R1(αp, φ_s, k3)
  R1p_turn = R1(αp, φ_o, k3)
  R1m_s = R1(αm, φ_s, k3)
  R1m_turn = R1(αm, φ_o, k3)

  Ip_s = -1/(B*rp2+A*rp1)*((B+A)*Fr_s + 2*r21*√(A*B)/(B*rp2-A*rp1)*R1p_s)
  Ip_turn = -1/(B*rp2+A*rp1)*((B+A)*Fr_turn + 2*r21*√(A*B)/(B*rp2-A*rp1)*R1p_turn)
  Im_s = -1/(B*rm2+A*rm1)*((B+A)*Fr_s + 2*r21*√(A*B)/(B*rm2-A*rm1)*R1m_s)
  Im_turn = -1/(B*rm2+A*rm1)*((B+A)*Fr_turn + 2*r21*√(A*B)/(B*rm2-A*rm1)*R1m_turn)

  Ip = Ip_turn - Ip_s
  Im = Im_turn - Im_s
  #if !isindir
    #Ip -= Ip_s
    #Im -= Im_s
  #else
  #  Ip += Ip_s
  #  Im += Im_s
  #end


  return real(2a/(rp-rm)*((rp-a*λ/2)*Ip - (rm-a*λ/2)*Im))
end

function f2(α, sinφ, j)
 p2 = √((1+α^2)/(1-j+α^2))
 return p2/2*log(abs((1-p2)/(1+p2)*(1+p2*√(1-j*sinφ^2))/(1-p2*√(1-j*sinφ^2))))
end
function S1(α, φ, j)
  return 1/(1+α^2)*(FastElliptic.F(φ, j) + α^2elliptic_pi(1+α^2, φ,j) - α*f2(α, sin(φ), j))
end
function I4ϕ(metric::Kerr, roots::NTuple{4}, root_diffs::NTuple{5}, rs, λ)
  a = metric.spin
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs
  rp = 1 + √(1-a^2)
  rm = 1 - √(1-a^2)

  if real(r32 * r41) < 0 || real(r31 * r42) < 0
    return 0
  end

  C = √real(r31 * r42)
  D = √real(r32 * r41)
  k4 = 4C * D / (C + D)^2
  a2 = abs(imag(r1))
  b1 = real(r4)

  k4 = 4 * C * D / (C + D)^2

  x4_s = (rs + b1) / a2
  x4_p = (rp + b1) / a2
  x4_m = (rm + b1) / a2


  go = √max((4a2^2 - (C - D)^2) / ((C + D)^2 - 4a2^2), 0.0)
  gp = (go*x4_p-1)/(go +x4_p)
  gm = (go*x4_m-1)/(go +x4_m)

  S1p_s = S1(gp, atan(x4_s) + atan(go), k4)
  S1p_full = S1(gp, π/2 + atan(go), k4)
  S1m_s = S1(gm, atan(x4_s) + atan(go), k4)
  S1m_full = S1(gm, π/2 + atan(go), k4)

  Fr_s = 2 / (C + D) * FastElliptic.F(atan(x4_s) + atan(go), k4)
  Fr_full = 2 / (C + D) * FastElliptic.F(π / 2 + atan(go), k4)

  Ip_s = go/(a2*(1-go*x4_p))*(Fr_s - 2/(C+D)*((1+go^2)/(go*(go+x4_p)))*S1p_s)
  Ip_full = go/(a2*(1-go*x4_p))*(Fr_full - 2/(C+D)*((1+go^2)/(go*(go+x4_p)))*S1p_full)
  Im_s = go/(a2*(1-go*x4_m))*(Fr_s - 2/(C+D)*((1+go^2)/(go*(go+x4_p)))*S1m_s)
  Im_full = go/(a2*(1-go*x4_m))*(Fr_full - 2/(C+D)*((1+go^2)/(go*(go+x4_p)))*S1m_full)

  Ip = Ip_full - Ip_s
  Im = Im_full - Im_s

  return real(2a/(rp-rm)*((rp-a*λ/2)*Ip - (rm-a*λ/2)*Im))
 
end

##----------------------------------------------------------------------------------------------------------------------
# θ Stuff
##----------------------------------------------------------------------------------------------------------------------
"""
  mino_time(η, λ, a, θs, θo, isindir::Bool, n::Int64)

Mino time of trajectory between two inclinations for a given screen coordinate

  `α` : Bardeen α

  `β` : Bardeen β 

  `a` : Blackhole angular Momentum

  `θs` : Emission inclination

  `θo` : Observer inclination

  `isindir` : Is the path direct or indirect?

  `n` : nth image in orde of amount of minotime traversed
"""
mino_time(metric::Kerr, α, β, θs, θo, isindir, n) = Gθ(metric, α, β, θs, θo, isindir, n)[1]

Gθ(metric::Kerr, α, β, θs, θo, isindir, n) = _Gθ(metric::Kerr, sign(β), θs, θo, isindir, n, η(metric, α, β, θo, metric.spin), λ(metric, α, θo))

function _Gθ(metric::Kerr, signβ, θs, θo, isindir, n, η, λ)
  a = metric.spin
  Go, Gs, Ghat, minotime, isvortical = 0.0, 0.0, 0.0, 0.0, η < 0.0

  isincone = abs(cos(θs)) < abs(cos(θo))
  if isincone && (isindir != ((signβ > 0) ⊻ (θo > π / 2)))
    return Inf, Inf, Inf, isvortical
  end
  if ((((signβ < 0) ⊻ (θs > π / 2)) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= π / 2) ⊻ (θs > π / 2)))
    return Inf, Inf, Inf, isvortical
  end

  Δθ = 1 / 2 * (1 - (η + λ^2) / a^2)
  up = Δθ + √(Δθ^2 + η / a^2)
  um = Δθ - √(Δθ^2 + η / a^2)
  m = up / um
  k = m

  #isvortical = η < 0.
  args = 0
  argo = 0
  k = 0
  if isvortical
    args = (cos(θs)^2 - um) / (up - um)
    argo = (cos(θo)^2 - um) / (up - um)
    k = 1.0 - m
    if (!(0.0 < argo < 1.0) || !(0.0 < args < 1.0))
      return Inf, Gs, Ghat, isvortical
    end
    tempfac = 1 / √abs(um * a^2)
    Go = ((θs > π / 2) ? -1 : 1) * tempfac * FastElliptic.F(asin(√argo), k)
    Gs = ((θs > π / 2) ? -1 : 1) * tempfac * FastElliptic.F(asin(√args), k)
    Ghat = 2tempfac * FastElliptic.K(k)
  else
    args = cos(θs) / √(up)
    argo = cos(θo) / √(up)
    k = m
    if !(-1 < args < 1) || !(-1 < argo < 1)
      return Inf, Gs, Ghat, isvortical
    end
    tempfac = 1 / √abs(um * a^2)
    #println("argo $(asin(argo)), $(k)")
    Go = tempfac * FastElliptic.F(asin(argo), k)
    #println("args $(asin(args)), $(k)")
    Gs = tempfac * FastElliptic.F(asin(args), k)
    Ghat = 2tempfac * FastElliptic.K(k)
  end

  νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > π / 2)
  minotime = real(isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
  return minotime, Gs, Ghat, isvortical
end

function Gs(metric::Kerr, α, β, θo, τ)
  return _Gs(metric, sign(β), θo, η(metric, α, β, θo, metric.spin), λ(metric, α, θo), τ)
end

function _Gs(metric::Kerr, signβ, θo, η, λ, τ)
  τ == Inf && return Inf
  a = metric.spin
  Go, Ghat, Ghat_2, isvortical = 0.0, 0.0, 0.0, η < 0.0

  Δθ = 1 / 2 * (1 - (η + λ^2) / a^2)
  up = Δθ + √(Δθ^2 + η / a^2)
  um = Δθ - √(Δθ^2 + η / a^2)
  m = up / um
  k = m

  argo = 0
  k = 0
  tempfac = 1 / √abs(um * a^2)

  if isvortical
    argo = (cos(θo)^2 - um) / (up - um)
    k = 1.0 - m
    tempfac = 1 / √abs(um * a^2)
    Go = tempfac * FastElliptic.F(asin(√argo), k)
    Ghat_2 = tempfac * FastElliptic.K(k)
    Ghat = 2Ghat_2
    Δτtemp = (τ%Ghat + (θo > π/2 ? -1 : 1)*signβ * Go)
    n = floor(τ/Ghat)
    Δτ = (θo > π/2 ? -1 : 1) * abs(argmin(abs, [(-1)^n*signβ*(Ghat - Δτtemp) , (-1)^n*signβ*Δτtemp]))
  else
    argo = cos(θo) / √(up)
    k = m
    tempfac = 1 / √abs(um * a^2)
    Go = tempfac * FastElliptic.F(asin(argo), k)
    Ghat_2 = tempfac * FastElliptic.K(k)
    Ghat = 2Ghat_2
    Δτtemp = (τ%Ghat + signβ * Go)
    n = floor(τ/Ghat)
    Δτ = argmin(abs, [(-1)^n*signβ*(Ghat - Δτtemp) , (-1)^n*signβ*Δτtemp])
  end

  return Δτ
end

function θs(metric::Kerr, α, β, θo, τ)
  return _θs(metric, sign(β), θo, η(metric, α, β, θo, metric.spin), λ(metric, α, θo), τ)
end

function _θs(metric::Kerr, signβ, θo, η, λ, τ)
  τ == Inf && return Inf, false
  a = metric.spin
  Go, Ghat, Ghat_2, isvortical = 0.0, 0.0, 0.0, η < 0.0

  Δθ = 1 / 2 * (1 - (η + λ^2) / a^2)
  up = Δθ + √(Δθ^2 + η / a^2)
  um = Δθ - √(Δθ^2 + η / a^2)
  m = up / um
  k = m

  #isvortical = η < 0.
  ans, k, argo = 0.0, 0.0, 0.0
  tempfac = 1 / √abs(um * a^2)
  if isvortical
    argo = (cos(θo)^2 - um) / (up - um)
    k = 1.0 - m
    tempfac = 1 / √abs(um * a^2)
    Go = tempfac * FastElliptic.F(asin(√argo), k)
    Ghat_2 = tempfac * FastElliptic.K(k)
    Ghat = 2Ghat_2
    Δτtemp = (τ%Ghat + (θo > π/2 ? -1 : 1)*signβ * Go)
    n = floor(τ/Ghat)
    absτs = abs(argmin(abs, [(-1)^n*signβ*(Ghat - Δτtemp) , (-1)^n*signβ*Δτtemp]))
    argr = (FastElliptic.sn(absτs/tempfac, k))^2
    ans = acos((θo > π/2 ? -1 : 1) * √((up-um)*argr + um))
  else
    argo = cos(θo) / √(up)
    k = m
    tempfac = 1 / √abs(um * a^2)
    Go = tempfac * FastElliptic.F(asin(argo), k)
    Ghat_2 = tempfac * FastElliptic.K(k)
    Ghat = 2Ghat_2
    Δτtemp = (τ%Ghat + signβ * Go)
    n = floor(τ/Ghat)
    τs = argmin(abs, [(-1)^n*signβ*(Ghat - Δτtemp) , (-1)^n*signβ*Δτtemp])
    #println(τs)
    newargs = FastElliptic.sn(τs/tempfac, k)
    ans = acos(√up*newargs)
  end

  return ans
end

Gϕ(metric::Kerr, α, β, θs, θo, isindir, n) = _Gϕ(metric::Kerr, sign(β), θs, θo, isindir, n, η(metric, α, β, θo, metric.spin), λ(metric, α, θo))

function _Gϕ(metric::Kerr, signβ, θs, θo, isindir, n, η, λ)
  a = metric.spin
  Go, Gs, Ghat, minotime, isvortical = 0.0, 0.0, 0.0, 0.0, η < 0.0

  isincone = abs(cos(θs)) < abs(cos(θo))
  if isincone && (isindir != ((signβ > 0) ⊻ (θo > π / 2)))
    return Inf, Inf, Inf, isvortical
  end
  if ((((signβ < 0) ⊻ (θs > π / 2)) ⊻ (n % 2 == 1)) && !isincone && !isvortical) || (isvortical && ((θo >= π / 2) ⊻ (θs > π / 2)))
    return Inf, Inf, Inf, isvortical
  end

  Δθ = 1 / 2 * (1 - (η + λ^2) / a^2)
  up = Δθ + √(Δθ^2 + η / a^2)
  um = Δθ - √(Δθ^2 + η / a^2)
  m = up / um
  k = m

  #isvortical = η < 0.
  args = 0
  argo = 0
  k = 0
  if isvortical
    args = (cos(θs)^2 - um) / (up - um)
    argo = (cos(θo)^2 - um) / (up - um)
    k = 1.0 - m
    if (!(0.0 < argo < 1.0) || !(0.0 < args < 1.0))
      return Inf, Gs, Ghat, isvortical
    end
    tempfac = 1 / ((1-um)*√abs(um * a^2))
    argn = (up-um)/(1-um)
    Go = ((θs > π / 2) ? -1 : 1) * tempfac * elliptic_pi(argn, asin(√argo), k)
    Gs = ((θs > π / 2) ? -1 : 1) * tempfac * elliptic_pi(argn, asin(√args), k)
    Ghat = 2tempfac * elliptic_pi(argn, k)
  else
    args = cos(θs) / √(up)
    argo = cos(θo) / √(up)
    k = abs(m)
    if !(-1 < args < 1) || !(-1 < argo < 1)
      return Inf, Gs, Ghat, isvortical
    end
    tempfac = 1 / √abs(um * a^2)
    Go = tempfac * elliptic_pi(up, asin(argo), k)
    Gs = tempfac * elliptic_pi(up, asin(args), k)
    Ghat = 2tempfac * elliptic_pi(up, k)
  end

  νθ = isincone ? (n % 2 == 1) ⊻ (θo > θs) : !isindir ⊻ (θs > π / 2)
  minotime = real(isindir ? (n + 1) * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs : n * Ghat - signβ * Go + (νθ ? 1 : -1) * Gs) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
  return minotime, Gs, Ghat, isvortical
end


##----------------------------------------------------------------------------------------------------------------------
# ϕ Stuff
##----------------------------------------------------------------------------------------------------------------------
function ϕs(metric::Kerr, α, β, θs, θo, rs, isindir, n)
  a = metric.spin
  ηtemp = η(metric, α, β, θo, a)
  λtemp = λ(metric, α, θo)

  roots = get_radial_roots(metric, ηtemp, λtemp)
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)

  root_diffs = get_root_diffs(roots...)
  Iϕ = 0
  if numreals == 4 #case2
    roots = real.(roots)
    root_diffs = real.(root_diffs)
    Iϕ = I2ϕ(metric, roots, root_diffs, rs, isindir, λtemp)
  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
      roots = (roots[1], roots[4], roots[2], roots[3])
    end
    Iϕ = I3ϕ(metric, roots, root_diffs, rs, isindir, λtemp)
  else
    Iϕ = I4ϕ(metric, roots, root_diffs, rs, λtemp)
  end
  Iϕ == Inf && return 0.0

  Gϕtemp,_,_,_ = Gϕ(metric, α, β, θs, θo, isindir, n)
  Gϕtemp == Inf && return 0.0
  return ϕs(Iϕ, Gϕtemp, λtemp, -π/2)
end

function ϕs(Iϕ, Gϕ, λ, ϕo)
  if Iϕ == NaN || Iϕ == Inf
    return 0.0
  end
  if Gϕ == NaN || Gϕ == Inf
    return 0.0
  end
  return Iϕ + λ*Gϕ +ϕo
  #return λ*Gϕ +ϕo
end

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
function p_boyer_lindquist_d(metric::Kerr, r, θ, η, λ, νr::Bool, νθ::Bool)
  a = metric.spin
  @SVector [-1, (νr ? 1 : -1) * √abs(r_potential(metric, η, λ, r)) / Δ(metric, r, a), λ, (νθ ? 1 : -1) * √abs(θ_potential(metric, η, λ, θ))]
end


"""
    jac_bl2zamo_du(r, θ, a)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
function jac_bl2zamo_du(metric::Kerr, r, θ)
  a = metric.spin
  # coords = {t, r, ϕ, θ}
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  At = A(metric, r, θ, a)

  return @SMatrix [# Eq 3.1 1972ApJ...178..347B
    √(At / (Σt * Δt)) 0.0         2*a*r/√(At * Σt * Δt) 0.0
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
    
    `a` : Blackhole spin
"""
function jac_zamo2bl_du(metric::Kerr, r, θ)
  a = metric.spin
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  At = A(metric, r, θ, a)

  return @SMatrix [
    # coords = {t, r, ϕ, θ}
    √((Σt * Δt) / At) 0.0                     -2*a*r*sin(θ)/√(At * Σt)  0.0
    0.0               √(Σt / Δt)              0.0                       0.0
    0.0               0.0                     √(At / Σt)*sin(θ)         0.0
    0.0               0.0                     0.0                       -√Σt
  ]
end

function jac_bl2zamo_ud(metric::Kerr, r, θ)
  a = metric.spin
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  At = A(metric, r, θ, a)

  return @SMatrix [#  Eq 3.2 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
    √((Σt * Δt) / At)             0.0         0.0               0.0
    0.0                           √(Σt / Δt)  0.0               0.0
    -(2a * r * sin(θ))/√(At * Σt) 0.0         √(At / Σt)*sin(θ) 0.0
    0.0                           0.0         0.0               -√Σt
  ]
end

function jac_zamo2bl_ud(metric::Kerr, r, θ)
  a = metric.spin
  Σt = Σ(metric, r, θ, a)
  Δt = Δ(metric, r, a)
  At = A(metric, r, θ, a)

  return @SMatrix [
    # coords = {t, r, ϕ, θ}
    √(At / (Σt * Δt))     0.0         0.0               0.0
    0.0                   √(Δt / Σt)  0.0               0.0
    2a*r/√(Σt * Δt * At)  0.0         √(Σt / At)*csc(θ) 0.0
    0.0                   0.0         0.0               -1/√Σt
  ]
end

function jac_zamo2fluid_ud(metric::Kerr, β, θ, φ)
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

function penrose_walker(metric::Kerr, r, θ, p_u::AbstractVector, f_u::AbstractVector)# Eq 6 arXiv:2001.08750v1
  a = metric.spin
  pt, pr, pϕ, pθ = p_u
  ft, fr, fϕ, fθ = f_u
  sinθ = sin(θ)
  cosθ = cos(θ)

  A = pt * fr - pr * ft + a * sinθ^2(pr * fϕ - pϕ * fr)
  B = ((r^2 + a^2) * (pϕ * fθ - pθ * fϕ) - a * (pt * fθ - pθ * ft)) * sinθ
  return A * r - B * a * cosθ, -(A * a * cosθ - B * r)
end

function screen_polarisation(metric::Kerr, κ::Complex, θ, α, β)# Eq 31 10.1103/PhysRevD.104.044060
  a = metric.spin
  #TODO: Check which is real and which is imaginary
  κ1 = real(κ)
  κ2 = imag(κ)

  μ = -(α + a * sin(θ))
  fα = (β * κ2 - μ * κ1) / (μ^2 + β^2)
  fβ = (β * κ1 + μ * κ2) / (μ^2 + β^2)


  return fα, fβ
end

evpa(fα, fβ) = atan(-fα, fβ)

function calcPol(metric::Kerr, α, β, ri, θs, θo, cross_spec_index, magfield::SVector{3, Float64}, βfluid::SVector{3, Float64}, νr::Bool, θsign::Bool)
  ri <= horizon(metric) && return 0.0, 0.0, 0.0, 0.0

  a = metric.spin
  βv = βfluid[1]
  θz = βfluid[2]
  ϕz = βfluid[3]

  ηtemp = η(metric, α, β, θo, a)
  λtemp = λ(metric, α, θo)
  p_bl_d = p_boyer_lindquist_d(metric, ri, θs, ηtemp, λtemp, νr, θsign)

  p_bl_u = metric_uu(metric, ri, θs) * p_bl_d
  p_zamo_u = jac_bl2zamo_ud(metric, ri, θs) * p_bl_u
  p_fluid_u = jac_zamo2fluid_ud(metric, βv, θz, ϕz) * p_zamo_u
  magfieldx, magfieldy, magfieldz = magfield
  _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
  vec = @SVector[p_fluid_uy * magfieldz - p_fluid_uz * magfieldy, p_fluid_uz * magfieldx - p_fluid_ux * magfieldz, p_fluid_ux * magfieldy - p_fluid_uy * magfieldx]
  norm = √sum(vec .* vec) + eps()
  f_fluid_u = SVector(zero(eltype(vec)), vec[1], vec[2], vec[3])
  f_zamo_u = jac_zamo2fluid_ud(metric, -βv, θz, ϕz) * f_fluid_u
  f_bl_u = jac_zamo2bl_ud(metric, ri, θs) * f_zamo_u
  A = @SMatrix [
    0.0   1.0           0.0         0.0
    -1.0  0.0           a*sin(θs)^2 0.0
    0.0   -a*sin(θs)^2  0.0         0.0
    0.0   0.0           0.0         0.0
  ]
  B = @SMatrix [
    0.0       0.0   0.0           -a*sin(θs)
    0.0       0.0   0.0           0.0
    0.0       0.0   0.0           (ri^2+a^2)
    a*sin(θs) 0.0   -(ri^2 + a^2) 0.0
  ]
  f_temp_d = ((A - B * im) * (ri - a * cos(θs) * im)) * (f_bl_u)
  κ = sum(p_bl_u .* f_temp_d)
  κ = κ / √(conj(κ) * κ)

  eα, eβ = screen_polarisation(metric, κ, θo, α, β) .* (norm^((cross_spec_index + 1.0) / 2))

  return eα, eβ, 1 / p_fluid_u[1], abs(p_fluid_u[1] / p_fluid_u[4])
end
