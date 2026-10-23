using Test
include(joinpath(@__DIR__, "gridRefinementStudy.jl"))

@testset "Unequal-ratio Richardson extrapolation and GCI" begin
    for hs in ((0.1, 0.2, 0.4), (0.1, 0.145, 0.22)), p in (0.5, 1.0, 1.7, 2.0), c in (-0.2, 0.2)
        values = 0.8 .+ c .* hs.^p
        fit = gci(hs..., values...)
        @test fit.status == "monotonic_fit"
        @test fit.p ≈ p atol=1e-10
        @test fit.extrapolated ≈ 0.8 atol=1e-12
        @test fit.uncertainty ≈ 1.25abs(c)*hs[1]^p rtol=1e-9
        @test fit.uncertainty_medium ≈ 1.25abs(c)*hs[2]^p rtol=1e-9
    end
    @test gci(1.0, 2.0, 4.0, 1.0, 0.9, 1.1).status == "oscillatory"
    @test gci(1.0, 2.0, 4.0, 1.0, 1.0, 1.0).status == "stagnant_or_roundoff"
    @test gci(1.0, 1.01, 2.0, 1.0, 1.1, 1.2).status == "ratios_too_small"
    @test gci(1.0, 2.0, 4.0, 1.0, 1.2, 1.3).status == "no_positive_order"
    @test isnan(gci(1.0, 2.0, 4.0, 0.0, 0.1, 0.3).relative_gci_pct)
end

@testset "Saved profile and field identity regression" begin
    root = normpath(joinpath(@__DIR__, ".."))
    paths = filter(p -> occursin(r"^verificationResults_\d+$", basename(p)),
                   readdir(root; join=true))
    for path in paths
        result = read_result(path)
        lines = readlines(joinpath(path, "log.verifyNormalShock"))
        for field in FIELDS
            index = findfirst(line -> strip(line)=="Field: $field", lines)
            @test !isnothing(index)
            # Restrict to this field's three norms: never read the next field.
            logged_l1 = metric(lines[index+1:index+3], "L1:")
            @test result.errors[(field, :point)][1] ≈ logged_l1 rtol=1e-7
        end
        @test result.h ≈ result.length_domain/result.n
        @test abs(result.xnum-result.xs) ≈ metric(lines, "Absolute error:") atol=1e-9
    end
end
