import Pkg
Pkg.activate((@__DIR__))
using Revise
using CairoMakie
using Jube

nmax = 1
spin = -0.8
θo = 10/180*π
θs = π/2
spec_index = 0.24
cross_spec_index = -1.5
rpeak = 4.7
p1 = 1.4
p2 = 5.7
βv = 0.97
χ = -2.19
ι = 1.3
o = Jube.AssymptoticObserver(1, θo)
acc = Jube.JuKeBOX(nmax, spin, spec_index, cross_spec_index, rpeak, p1, p2, βv, χ, ι, χ+π)
function raytrace_point(acc, α, β, θs, o)
    direct  = Jube.raytrace(acc, -α, β, θs, o, true)
    indirect = Jube.raytrace(acc, -α, β, θs, o, false)
    return direct .+ indirect
end

sze = 200
αvals = range(-20, 20, length=sze)
βvals = range(-20, 20, length=sze)
ivals = zeros(sze, sze)
qvals = zeros(sze, sze)
uvals = zeros(sze, sze)
vvals = zeros(sze, sze)

for i in 1:sze
    for j in 1:sze
        stokes = raytrace_point(acc, αvals[i], βvals[j], θs, o)
        ivals[i, j] = stokes[1]
        qvals[i, j] = stokes[2]
        uvals[i, j] = stokes[3]
        vvals[i, j] = stokes[4]

    end
end

fig = Figure();
ax0 = Axis(fig[1, 1], aspect=1, title="Stokes I")
ax1 = Axis(fig[1, 3], aspect=1, title="Stokes Q")
ax2 = Axis(fig[2, 1], aspect=1, title="Stokes U")
ax3 = Axis(fig[2, 3], aspect=1, title="Stokes V")

hm0 = CairoMakie.heatmap!(ax0, ivals, colormap=:afmhot)
Colorbar(fig[1, 2], hm0)
hm1 = CairoMakie.heatmap!(ax1, qvals, colormap=:afmhot)
Colorbar(fig[1, 4], hm1)
hm2 = CairoMakie.heatmap!(ax2, uvals, colormap=:afmhot)
Colorbar(fig[2, 2], hm2)
hm3 = CairoMakie.heatmap!(ax3, vvals, colormap=:afmhot)
Colorbar(fig[2, 4], hm3)


display(fig)
