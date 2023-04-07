function _make_pol_data_field(met::Kerr, αmax, βmax, n, steps, θs, o, B, βfluid)

    rh = horizon(met)
    θo = o.inclination


    αvals = LinRange(-αmax,αmax, 2steps)
    βvals = LinRange(-βmax, βmax, 2steps)
    rsvals = [zeros(2steps) for _ in 1:2steps]
    rsvals2 = [zeros(2steps) for _ in 1:2steps]
    rootvals = [zeros(2steps) for _ in 1:2steps]
    sinvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
    cosvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]

    isincone = abs(cos(θs)) < abs(cos(θo))
    νθtrue =  isincone ? (n%2==1) ⊻ (θo > θs) : false ⊻ (θs > π/2)
    νθfalse =  isincone ? (n%2==1) ⊻ (θo > θs) : true ⊻ (θs > π/2)
    
    for _i in 1:(2steps)
        α = αvals[_i]
        for _j in 1:(2steps)
            β = βvals[_j]
            rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, true, n)
            evalpol = (rstemp >= rh && rstemp != Inf)
            sintemp, costemp =  evalpol ? Jube.calcPol(met, α, β, rstemp, θs, θo, 1.0, B, βfluid, νrtemp, νθtrue) : (0., 0.)

            rsvals[_i][_j] = rstemp
            rootvals[_i][_j] = numroots
            sinvals[_i][_j] = sintemp * αmax 
            cosvals[_i][_j] = costemp * βmax

            rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, false, n)
            evalpol = (rstemp >= rh && rstemp != Inf)
            sintemp, costemp =  evalpol ? Jube.calcPol(α, β, rstemp, θs, θo, 1.0, B, βfluid, νrtemp, νθfalse) : (0., 0.)

            rsvals2[_i][_j] += rstemp
            rootvals[_i][_j] = numroots
            sinvals[_i][_j] = sintemp * αmax
            cosvals[_i][_j] = costemp * βmax
        end
    end
    if abs(cos(θo)) > abs(cos(θs))
        rsvals .+= rsvals2
    end

    return rsvals, rootvals, sinvals, cosvals
end

function _make_radial_data_field(met::Kerr, αmax, βmax, n, steps, θs, o) 
    θo = o.inclination

    αvals = LinRange(-αmax,αmax, 2steps)
    βvals = LinRange(-βmax, βmax, 2steps)
    rsvals = [zeros(2steps) for _ in 1:2steps]
    rsvals2 = [zeros(2steps) for _ in 1:2steps]
    rootvals = [zeros(2steps) for _ in 1:2steps]
    
    for _i in 1:(2steps)
        α = αvals[_i]
        for _j in 1:(2steps)
            β = βvals[_j]
            rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, true, n)

            rsvals[_i][_j] = rstemp
            rootvals[_i][_j] = numroots

            rstemp, νrtemp, numroots = Jube.rs(met, α, β, θs, o, false, n)

            rsvals2[_i][_j] += rstemp
            rootvals[_i][_j] = numroots
        end
    end
    if abs(cos(θo)) > abs(cos(θs))
        rsvals .+= rsvals2
    end

    return rsvals, rootvals
end

# Abstract Metric Interface Tests
@testset "abstract metric interface" begin
    #test metric interface
    struct TestMetric <: Jube.AbstractMetric end
    @test_throws MethodError Jube.metric_dd(TestMetric())
    @test_throws MethodError Jube.metric_uu(TestMetric())

    Jube.metric_dd(met::TestMetric, r) = [1 0; 0 r]
    @test Jube.metric_uu(TestMetric(), 10) == [1 0; 0 1/10]
end


# Kerr tests
@testset "horizon" begin
    @test horizon(Jube.Kerr(0.0)) == (2)
    @test horizon(Jube.Kerr(0.5)) == (1 + √(1 - 0.5^2))
    @test horizon(Jube.Kerr(1.0)) == (1)
end

@testset "get_roots" begin
    # test Schwarzschild
    @test all(isapprox.(Jube.get_radial_roots(Kerr(0.0), 27, 0.), (-6, 0, 3, 3), rtol=1e-5))
    @test all(isapprox.(Jube.get_radial_roots(Kerr(0.0), 0., √27.), (-6, 0, 3, 3), rtol=1e-5))
    # test case 1 or 2
    @test all(isapprox.(Jube.get_radial_roots(Kerr(-0.82), -1.20503, 11.4), (-12.332812618682777, -0.002731857389424548, 2.4271811157108694, 9.908363360361331), rtol=1e-5))
    # test case 3
    @test begin 
        roots = Jube.get_radial_roots(Kerr(-0.915), 2.47908, 2.4)
        roots_ans = (-3.820331170598019, 0.07875611823150777, 1.8707875261832558 - 1.8435163296601766im, 1.8707875261832558 + 1.8435163296601766im )
        (roots[3] ≈ conj(roots[4])) && all(isapprox.(map(x->real(x)+abs(imag(x))im, roots), map(x->real(x)+abs(imag(x))im, roots_ans), rtol=1e-5))
        
    end
    # test for case 4
    @test begin
        roots = Jube.get_radial_roots(Kerr(-1.0),-0.087975, -1)
        roots_ans = (-0.3710972685097377 - 0.5479363974641974im, -0.3710972685097377 + 0.5479363974641974im, 0.3710972685097377 - 0.2513305984649188im, 0.3710972685097377 + 0.2513305984649188im)
        (roots[1] ≈ conj(roots[2])) && (roots[3] ≈ conj(roots[4]))all(isapprox.(map(x->real(x)+abs(imag(x))im, roots), map(x->real(x)+abs(imag(x))im, roots_ans), rtol=1e-5))
    end
end

@testset "emission coordinates" begin
    ##Test θs inversion
    #@test begin 
    #    check = true
    #    for spin in [-0.95, 0.5, 0.01]
    #        if check == false
    #            break
    #        end
    #        met = Jube.Kerr(spin)
    #        all_vals = []

    #        for j in 1:1
    #            θo = j*π/8
    #            for i in 1:1:179
    #                tempθs = i/180*π
    #                if θo*180/π == i || θo*180/π == (180 - i) || i == 90
    #                    continue
    #                end
    #                for isindir in [true, false]
    #                    for n in [0, 1, 2]
    #                        for β in -10.0:0.1:10.0
    #                            #TODO: Fix behaviour of raytracing for β = 0
    #                            if isapprox(β, 0.0, atol=1e-6)
    #                                continue
    #                            end
    #                            for α in -10.0:0.1:10.0
    #                                τ, τs, _, _ = Jube.Gθ(met, α, β, tempθs, θo, isindir, n)
    #                                temp = Jube.Gs(met, α, β, θo, τ) 
    #                                if τ == Inf
    #                                    append!(all_vals,true)
    #                                else
    #                                    sol = isapprox(temp, τs, rtol=1e-6)
    #                                    append!(all_vals, sol)
    #                                end
    #                            end
    #                        end
    #                    end
    #                end
    #            end
    #        end
    #        check = all(all_vals)
    #    end
    #    check
    #end
    @test begin 
        check = true
        for spin in [-0.99, 0.5, 0.01]
            if check == false
                break
            end
            met = Jube.Kerr(spin)
            all_vals = []

            for j in 1:1
                θo = j*π/8
                for i in 1:1:179
                    tempθs = i/180*π
                    if θo*180/π == i || θo*180/π == (180 - i) || i == 90
                        continue
                    end
                    for isindir in [true, false]
                        for n in [0, 1, 2]
                            for β in -10.0:0.1:10.0
                                #TODO: Fix behaviour of raytracing for β = 0
                                if isapprox(β, 0.0, atol=1e-6)
                                    continue
                                end
                                for α in -10.0:0.1:10.0
                                    τ, τs, _, _ = Jube.Gθ(met, α, β, tempθs, θo, isindir, n)
                                    temp = Jube.θs(met, α, β, θo, τ) 
                                    if τ == Inf
                                        append!(all_vals,true)
                                    else
                                        sol = isapprox(temp, tempθs, rtol=1e-6)
                                        append!(all_vals, sol)
                                    end
                                end
                            end
                        end
                    end
                end
            end
            check = all(all_vals)
        end
        check
    end
    #@test begin 
    #    check = true
    #    for spin in [-0.95, 0.5, 0.01]
    #        if check == false
    #            break
    #        end
    #        met = Jube.Kerr(spin)
    #        all_vals = []

    #        for j in 1:1
    #            θo = j*π/8
    #            o = Jube.AssymptoticObserver(1, θo)
    #            for i in 1:1:179
    #                tempθs = i/180*π
    #                if θo*180/π == i || θo*180/π == (180 - i) || i == 90
    #                    continue
    #                end
    #                for isindir in [true, false]
    #                    for n in [0, 1, 2]
    #                        for β in -10.0:0.1:10.0
    #                            #TODO: Fix behaviour of raytracing for β = 0
    #                            if isapprox(β, 0.0, atol=1e-6)
    #                                continue
    #                            end
    #                            for α in -10.0:0.1:10.0
    #                                ηtemp = Jube.η(met, α, β, θo, met.spin)
    #                                λtemp = Jube.λ(met, α, θo)

    #                                τ, τs, _, _ = Jube.Gθ(met, α, β, tempθs, θo, isindir, n)
    #                                rs, νr, numroots = Jube._rs(met, ηtemp, λtemp, τ)
    #                                
    #                                if τ == Inf || rs < horizon(met)
    #                                    append!(all_vals,true)
    #                                else
    #                                    τr = Jube.Ir(met, νr, θo ,rs, α, β)
    #                                    sol = isapprox(τr, τ, rtol=1e-6)
    #                                    append!(all_vals, sol)
    #                                end
    #                            end
    #                        end
    #                    end
    #                end
    #            end
    #        end
    #        check = all(all_vals)
    #    end
    #    check
    #end

end
        
#    for radialTestData in readdir((@__DIR__)*"/radialTestData")
#        path = (@__DIR__)*"/radialTestData/"*radialTestData
#        data = JSON.parsefile(path)
#        αmax = data["αmax"]
#        βmax = data["βmax"]
#        steps = data["steps"]
#        rsvals = data["rsvals"]
#        θo = data["θo"]
#        θs = data["θs"]
#        a = data["a"]
#        n = data["n"]
#        rsvals = data["rsvals"]
#        rootsvals = data["rootsvals"]
#
#        met = Kerr(a)
#        o = Jube.AssymptoticObserver(1, θo)
#        testrsvals, testrootsvals = _make_radial_data_field(met, αmax, βmax, n, steps, θs, o)
#        @test begin
#             (abs(maximum(maximum.(rsvals .- testrsvals))) < eps())
#        end
#        @test begin
#             (abs(maximum(maximum.(rootsvals .- testrootsvals))) < 1)
#        end
#    end
#end

# @testset "emission_pol_field" begin
#    for polTestData in readdir((@__DIR__)*"/polTestData")
#        path = (@__DIR__)*"/polTestData/"*polTestData
#        data = JSON.parsefile(path)
#        αmax = data["αmax"]
#        βmax = data["βmax"]
#        steps = data["steps"]
#        θo = data["θo"]
#        θs = data["θs"]
#        a = data["a"]
#        n = data["n"]
#        βfluid = data["βfluid"]
#        B = data["B"]
#        sinvals = data["sinvals"]
#        cosvals = data["cosvals"]
#
#        met = Kerr(a)
#        o = Jube.AssymptoticObserver(1, θo)
#        _, _, testsinvals, testcosvals = _make_pol_data_field(met, αmax, βmax, n, steps, θs, o, [i for i in B], [i for i in βfluid])
#
#        @test begin
#             (abs(maximum(maximum.(sinvals .- testsinvals))) < eps())
#        end
#        @test begin
#             (abs(maximum(maximum.(cosvals .- testcosvals))) < eps())
#        end
#
#    end
#end