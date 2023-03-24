import Pkg
Pkg.activate((@__DIR__))
using Revise
using CairoMakie
using Jube
using StaticArrays
#include((@__DIR__)*"/src/Jube.jl")

function plot(αmax, βmax, rmin, rmax, ntemp, a, steps, θs, θo, B, βfluid)
    met = Jube.Kerr(a)
    o = Jube.AssymptoticObserver(1, θo)

    rh = 1 + √(1-a^2)
    
    f = Figure(resolution = (800, 800))
    ax0 = Axis(f[1,1], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    #xlims!(ax0, -αmax, αmax)
    #ylims!(ax0, -βmax, βmax)

    ax1 = Axis(f[1,2], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    #xlims!(ax1, -αmax, αmax)
    #ylims!(ax1, -βmax, βmax)

    ax2 = Axis(f[2,1], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    #xlims!(ax2, -αmax, αmax)
    #ylims!(ax2, -βmax, βmax)
    
    ax3 = Axis(f[2,2], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    #xlims!(ax3, -αmax, αmax)
    #ylims!(ax3, -βmax, βmax)
    αvalstemp = range(-αmax,αmax, length=2steps)
    βvalstemp = range(-βmax, βmax, length=2steps)
    sinvals =     [zeros(2steps) for _ in 1:2steps]
    cosvals =     [zeros(2steps) for _ in 1:2steps]
    for n in 0:ntemp


        αvals = (cos(pa) .* αvalstemp) .+ (sin(pa) .* βvalstemp)
        βvals = (cos(pa) .* βvalstemp) .- (sin(pa) .* αvalstemp)

        rsvals =      [zeros(2steps) for _ in 1:2steps]
        rsvals2 =     [zeros(2steps) for _ in 1:2steps]
        rsvals3 =     [zeros(2steps) for _ in 1:2steps]
        rsvals4 =     [zeros(2steps) for _ in 1:2steps]

        geotypevals = [zeros(2steps) for _ in 1:2steps]
        rootvals =    [zeros(2steps) for _ in 1:2steps]

        #νθ = θo > θs 

        isincone = abs(cos(θs)) < abs(cos(θo))
        νθtrue =  isincone ? (n%2==1) ⊻ (θo > θs) : false ⊻ (θs > π/2)
        νθfalse =  isincone ? (n%2==1) ⊻ (θo > θs) : true ⊻ (θs > π/2)
        for _i in 1:(2steps)
            α = αvalstemp[_i]
            for _j in 1:(2steps)
                β = βvalstemp[_j]
                rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, true, n)
                evalpol = (rstemp >= rh && rstemp != Inf)
                sintemp, costemp =  evalpol ? Jube.calcPol(met, α, β, rstemp, θs, θo, 1, B, βfluid, νrtemp, νθtrue) : (0., 0.)

                rstemp2, νrtemp2, numroots2 = Jube.rs(met, α, β, π-θs, o, true, n)
                evalpol = (rstemp2 >= rh && rstemp2 != Inf)
                sintemp, costemp =  evalpol ? Jube.calcPol(met, α, β, rstemp, π-θs, θo, 1, B, βfluid, νrtemp, νθtrue) : (0., 0.)

                rsvals[_i][_j] = rstemp
                rsvals3[_i][_j] = rstemp2
                rootvals[_i][_j] = numroots
                sinvals[_i][_j] += sintemp * αmax / 20
                cosvals[_i][_j] += costemp * βmax / 20

                rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, false, n)
                evalpol = (rstemp >= rh && rstemp != Inf)
                sintemp, costemp =  evalpol ? Jube.calcPol(met, α, β, rstemp, θs, θo, 1, B, βfluid, νrtemp, νθtrue) : (0., 0.)

                rstemp2, νrtemp2, numroots2 = Jube.rs(met, α, β, π-θs, o, false, n)
                evalpol = (rstemp2 >= rh && rstemp2 != Inf)
                sintemp, costemp =  evalpol ? Jube.calcPol(met, α, β, rstemp, π-θs, θo, 1, B, βfluid, νrtemp, νθtrue) : (0., 0.)

                rsvals2[_i][_j] += rstemp
                rsvals4[_i][_j] += rstemp2
                rootvals[_i][_j] = numroots
                sinvals[_i][_j] += sintemp * αmax / 20
                cosvals[_i][_j] += costemp * βmax / 20

                geotypevals[_i][_j] = (Jube.Gθ(met, α, β, θs, θo, false, n)[3] ? 1 : 0)
            end
        end
        tickx = reduce(vcat,transpose.(sinvals)) 
        ticky = reduce(vcat,transpose.(cosvals))

        arrows!(ax2, .-αvalstemp, βvalstemp, .-tickx, ticky, arrowhead = ' ')

        if abs(cos(θo)) > abs(cos(θs))
            rsvals .+= rsvals2
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/5:rmax)
        else 
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/5:rmax)
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals2)), levels=rmin:(rmax - rmin)/5:rmax)
        end

        if abs(cos(θo)) > abs(cos(θs))
            rsvals3 .+= rsvals4
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals3)),  levels=rmin:(rmax - rmin)/5:rmax)
        else 
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals3)),  levels=rmin:(rmax - rmin)/5:rmax)
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals4)), levels=rmin:(rmax - rmin)/5:rmax)
        end

        heatmap!(ax0, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rootvals)))
        heatmap!(ax3, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(geotypevals)))

    end
    display(f)
end

αmax = 20
βmax = 20

a = 0.92872191 
θo = 0.45379006
pa = 0.60225741
spec_index = 0.82056789
rpeak = 0.28258731
p1 = 0.49620015
p2 = 0.87419257
βv = 0.8807092
ϕz = 0.6482512
ι = 0.47341764
steps = 50
n = 2
θs = (90)*π/180
back = true

θz = π/2
B = @SVector[sin(ι)*cos(ϕz), sin(ι)*sin(ϕz), cos(ι)]
βfluid = @SVector[βv, θz, ϕz]
rmin = 2
rpeak = 3.5
rmax = 10

count = 0
@time for i in range( 1, 179, length=179)
    θo = i*π/180
    if θo == θs || (180-i)*π/180 == θs
        continue
    end
    plot(αmax, βmax, rmin, rmax, n, a, steps, θs, θo, B, βfluid)

    count +=1
end
