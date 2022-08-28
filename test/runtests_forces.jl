
@testset "Forces" begin
    # From potential

    force = ForceFromPotential("Flat")
    @test force_eval(force, [0.0])[1] ≈ 0.0
    @test force_eval(force, [1.0])[1] ≈ 0.0

    force = ForceFromPotential("Flat",3)
    @test force_eval(force, [0.0,0.0,0.0]) ≈ [0.0,0.0,0.0]
    @test force_eval(force, [1.0,0.0,0.0]) ≈ [0.0,0.0,0.0]

    force = ForceFromPotential("Harmonic")
    @test force_eval(force, [0.0])[1] ≈ 0.0
    @test force_eval(force, [1.0])[1] ≈ -1.0

    force = ForceFromPotential("DoubleWell")
    @test force_eval(force, [0.0])[1] ≈ 0.0
    @test force_eval(force, [1.0])[1] ≈ 0.0

    force=ForceFromPotential("Muller", 2)
    @test force_eval(force,[0.0,0.0]) ≈ [-120.44528523713869, -108.79148986312214]
    @test force_eval(force,[-0.70,1.25]) ≈ [ 32.34319126842259, -146.34307531443176]

    #Forces from ApproxFun
    force = ForceFromBasis("Taylor", [0.0 -2.0])
    @test force_eval(force, [0.0])[1] ≈ 0.0
    @test force_eval(force, [1.0])[1] ≈ -2.0

    #Forces from BSplinesKit
    force = ForceFromSplines(
        3,
        [0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],
        [8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538],
    )

    @test force_eval(force, [0.5])[1] ≈ -0.07372447
    @test force_eval(force, [1.0])[1] ≈ -9.0294072
    # @test force_eval(force,[2.0])[1] ≈ 49.3164884 # BSplineKit does not extrapolate TODO

    #Forces from Scipy.interpolate
    force = ForceFromScipySplines(
        3,
        [0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],
        [8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538],
    )

    @test force_eval(force, [0.5])[1] ≈ -0.07372447
    @test force_eval(force, [1.0])[1] ≈ -9.0294072
    @test force_eval(force, [2.0])[1] ≈ 49.3164884
end
