using Printf
using Statistics
using DelimitedFiles

const FIELDS = (:p, :rho, :T, :u)
const PROFILE_HEADER = "x,p_cfd,p_exact,rho_cfd,rho_exact,T_cfd,T_exact,u_cfd,u_exact"

function options(args)
    case_dir = normpath(joinpath(@__DIR__, ".."))
    field_tol, shock_tol, gci_tol = 1e-3, 1e-4, 1e-3
    for arg in args
        if startswith(arg, "--field-tol=")
            field_tol = parse(Float64, split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--shock-tol=")
            shock_tol = parse(Float64, split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--gci-tol=")
            gci_tol = parse(Float64, split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            case_dir = abspath(arg)
        end
    end
    all(t -> isfinite(t) && t > 0, (field_tol, shock_tol, gci_tol)) ||
        error("Tolerances must be finite positive fractions, not percentages")
    return (; case_dir, field_tol, shock_tol, gci_tol)
end

function metric(lines, label)
    matches = filter(line -> startswith(strip(line), label), lines)
    length(matches) == 1 || error("Expected one '$label' entry; found $(length(matches))")
    return parse(Float64, first(split(strip(split(only(matches), label; limit=2)[2]))))
end

function crossing(x, q, target, near)
    candidates = Float64[]
    for i in 1:length(x)-1
        if q[i] != q[i+1] && (q[i] - target) * (q[i+1] - target) <= 0
            push!(candidates, x[i] + (target-q[i]) * (x[i+1]-x[i]) / (q[i+1]-q[i]))
        end
    end
    isempty(candidates) && error("No pressure crossing at $target")
    return candidates[argmin(abs.(candidates .- near))]
end

norms(e) = (mean(abs.(e)), sqrt(mean(abs2, e)), maximum(abs.(e)))
rate(ec, ef, hc, hf) = ec > 0 && ef > 0 ? log(ec/ef)/log(hc/hf) : NaN

function read_result(dir)
    n = parse(Int, match(r"^verificationResults_(\d+)$", basename(dir)).captures[1])
    lines = readlines(joinpath(dir, "log.verifyNormalShock"))
    metric(lines, "Cells:") == n || error("Log/directory cell count mismatch: $dir")
    time = metric(lines, "Time:")
    xs = metric(lines, "Exact shock location:")
    path = joinpath(dir, "verification_profiles.csv")
    first(readlines(path)) == PROFILE_HEADER || error("Unexpected profile header: $path")
    data = readdlm(path, ',', Float64; skipstart=1)
    size(data) == (n, 9) || error("Profile size/cell count mismatch: $path")
    all(isfinite, data) || error("Nonfinite profile values: $path")
    data = data[sortperm(data[:, 1]), :]
    x = data[:, 1]
    h = mean(diff(x))
    h > 0 && maximum(abs.(diff(x) .- h)) < 1e-6*h ||
        error("Only uniform 1D grids are supported: $dir")
    left, right = x[1]-h/2, x[end]+h/2
    length_domain = n*h
    left < xs < right || error("Shock is outside the domain: $dir")
    p = data[:, 2]
    p2, p1 = data[1, 3], data[end, 3]
    jump_p = p2-p1
    xnum = crossing(x, p, p1+0.5jump_p, xs)
    thickness = abs(crossing(x, p, p1+0.9jump_p, xs) -
                    crossing(x, p, p1+0.1jump_p, xs))
    # Signed error removes the exact motion at each run's actual output time.
    qois = Dict("shock_offset" => (xnum-xs)/length_domain)
    scales = Dict("shock_offset" => length_domain)
    errors = Dict{Tuple{Symbol, Symbol}, NTuple{3, Float64}}()
    mask2, mask1 = x .< xs-0.03, x .> xs+0.03
    any(mask1) && any(mask2) || error("No plateau samples: $dir")
    fraction = clamp.((xs .- (x .- h/2)) ./ h, 0, 1)
    for (j, field) in enumerate(FIELDS)
        numerical, exact_point = data[:, 2j], data[:, 2j+1]
        q2, q1 = exact_point[1], exact_point[end]
        jump = abs(q2-q1)
        jump > 0 || error("Zero jump for $field: $dir")
        exact_average = q1 .+ (q2-q1) .* fraction
        errors[(field, :point)] = norms(numerical .- exact_point)
        errors[(field, :cell_average)] = norms(numerical .- exact_average)
        errors[(field, :normalized_point)] = norms((numerical .- exact_point)./jump)
        errors[(field, :normalized_cell_average)] = norms((numerical .- exact_average)./jump)
        # Domain-mean defect and post-shock plateau defect, scaled by the jump.
        qois["$(field)_mean_defect"] = mean(numerical .- exact_average)/jump
        qois["$(field)_plateau_defect"] = (mean(numerical[mask2])-q2)/jump
        scales["$(field)_mean_defect"] = jump
        scales["$(field)_plateau_defect"] = jump
    end
    solver_log = joinpath(dir, "log.rhoCentralFoam")
    runtime = NaN
    if isfile(solver_log)
        executions = collect(eachmatch(r"ExecutionTime = ([0-9.eE+-]+)", read(solver_log, String)))
        !isempty(executions) && (runtime = parse(Float64, last(executions).captures[1]))
    end
    return (; n, h, time, xs, xnum, thickness, length_domain, left, runtime,
             qois, scales, errors)
end

# Solve the unequal-ratio power-law model phi(h)=phi0+C*h^p.
# Fine=1, medium=2, coarse=3. Oscillating data do not support this model.
function gci(h1, h2, h3, f1, f2, f3)
    r21, r32 = h2/h1, h3/h2
    d21, d32 = f2-f1, f3-f2
    empty_result(status) = (; status, p=NaN, extrapolated=NaN, uncertainty=NaN,
                            uncertainty_medium=NaN, relative_gci_pct=NaN)
    min(r21, r32) >= 1.1 || return empty_result("ratios_too_small")
    min(abs(d21), abs(d32)) > 1e-12*max(1, abs(f1), abs(f2), abs(f3)) ||
        return empty_result("stagnant_or_roundoff")
    sign(d21) == sign(d32) || return empty_result("oscillatory")
    residual(p) = p*log(r21) + log(expm1(p*log(r32))) -
                  log(expm1(p*log(r21))) - log(abs(d32/d21))
    lo, hi = 1e-4, 12.0
    residual(lo) < 0 < residual(hi) || return empty_result("no_positive_order")
    for _ in 1:100
        mid = (lo+hi)/2
        if residual(mid) > 0
            hi = mid
        else
            lo = mid
        end
    end
    p = (lo+hi)/2
    correction = (f1-f2)/expm1(p*log(r21))
    uncertainty = 1.25abs(correction)
    uncertainty_medium = 1.25abs(d32)/expm1(p*log(r32))
    relative_gci_pct = abs(f1) > 1e-12 ? 100uncertainty/abs(f1) : NaN
    return (; status="monotonic_fit", p, extrapolated=f1+correction,
            uncertainty, uncertainty_medium, relative_gci_pct)
end

function write_table(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(keys(first(rows)), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function study(config)
    names = filter(name -> occursin(r"^verificationResults_\d+$", name) &&
                   isdir(joinpath(config.case_dir, name)), readdir(config.case_dir))
    length(names) >= 3 || error("At least three saved grid levels are required")
    grids = sort(read_result.(joinpath.(config.case_dir, names)); by=g -> g.n)
    all(g -> isapprox(g.length_domain, grids[1].length_domain; rtol=1e-6) &&
             isapprox(g.left, grids[1].left; atol=1e-8), grids) || error("Domains differ")
    for name in keys(grids[1].scales)
        all(g -> isapprox(g.scales[name], grids[1].scales[name]; rtol=1e-6), grids) ||
            error("Reference states differ between grids: $name")
    end
    out = joinpath(config.case_dir, "gridRefinementResults")
    mkpath(out)
    mesh_rows, error_rows, gci_rows = NamedTuple[], NamedTuple[], NamedTuple[]
    for (i, g) in enumerate(grids)
        push!(mesh_rows, (; cells=g.n, dx_m=g.h, time_s=g.time,
            refinement_ratio=i==1 ? NaN : grids[i-1].h/g.h,
            shock_error_m=abs(g.xnum-g.xs), shock_error_cells=abs(g.xnum-g.xs)/g.h,
            shock_thickness_m=g.thickness, shock_thickness_cells=g.thickness/g.h,
            solver_seconds=g.runtime))
        for field in FIELDS, basis in (:point, :cell_average), k in 1:3
            e = g.errors[(field, basis)][k]
            en = g.errors[(field, Symbol("normalized_", basis))][k]
            order = i==1 ? NaN : rate(grids[i-1].errors[(field, basis)][k], e,
                                      grids[i-1].h, g.h)
            push!(error_rows, (; cells=g.n, field, basis, norm=("L1", "L2", "Linf")[k],
                               error=e, error_over_jump=en, pairwise_order=order))
        end
    end
    # Stride 2 provides ~1.44 ratios in addition to the ~1.2 adjacent ratios.
    for stride in (1, 2), i in 2stride+1:length(grids), name in sort(collect(keys(grids[i].qois)))
        fine, medium, coarse = grids[i], grids[i-stride], grids[i-2stride]
        f1, f2, f3 = fine.qois[name], medium.qois[name], coarse.qois[name]
        fit = gci(fine.h, medium.h, coarse.h, f1, f2, f3)
        check_ratio = NaN
        # A fourth grid tests the fitted power law independently of its three inputs.
        if i > 3stride && isfinite(fit.p)
            fourth = grids[i-3stride]
            predicted = (f3-f2) * ((fourth.h/coarse.h)^fit.p-1) /
                        (1-(medium.h/coarse.h)^fit.p)
            check_ratio = (fourth.qois[name]-f3)/predicted
        end
        supported = fit.status == "monotonic_fit" && fit.p <= 2 &&
                    isfinite(check_ratio) && 0.8 <= check_ratio <= 1.25
        covered = isfinite(fit.uncertainty) && abs(f1) <= fit.uncertainty
        scale = fine.scales[name]
        push!(gci_rows, (; quantity=name, stride, fine_cells=fine.n,
            medium_cells=medium.n, coarse_cells=coarse.n,
            r21=medium.h/fine.h, r32=coarse.h/medium.h,
            fine_value_scaled=f1, exact_value_scaled=0.0, observed_order=fit.p,
            extrapolated_scaled=fit.extrapolated,
            fine_gci_scaled=fit.uncertainty, medium_gci_scaled=fit.uncertainty_medium,
            fine_gci_physical=fit.uncertainty*scale,
            scale_physical=scale, relative_gci_percent=fit.relative_gci_pct,
            lower_scaled=f1-fit.uncertainty, upper_scaled=f1+fit.uncertainty,
            fourth_grid_ratio=check_ratio, exact_covered=covered,
            supported, status=fit.status))
    end
    assessments = NamedTuple[]
    for g in grids
        max_l1 = maximum(g.errors[(f, :normalized_cell_average)][1] for f in FIELDS)
        shock_error = abs(g.qois["shock_offset"])
        rows = filter(r -> r.fine_cells == g.n && r.stride == 1, gci_rows)
        supported = length(rows)==length(g.qois) && all(r -> r.supported && r.exact_covered, rows)
        max_gci = isempty(rows) || !all(r -> isfinite(r.fine_gci_scaled), rows) ?
                  NaN : maximum(r.fine_gci_scaled for r in rows)
        accuracy_pass = max_l1 <= config.field_tol && shock_error <= config.shock_tol
        push!(assessments, (; cells=g.n, dx_m=g.h, max_field_L1_over_jump=max_l1,
            shock_error_over_length=shock_error, max_gci_scaled=max_gci,
            accuracy_pass, gci_supported=supported,
            eligible=accuracy_pass && supported && max_gci <= config.gci_tol))
    end
    write_table(joinpath(out, "grid_refinement_summary.csv"), mesh_rows)
    write_table(joinpath(out, "field_errors.csv"), error_rows)
    write_table(joinpath(out, "gci_triplets.csv"), gci_rows)
    write_table(joinpath(out, "mesh_selection.csv"), assessments)
    report(config, grids, gci_rows, assessments, out)
    return (; grids, error_rows, gci_rows, assessments, out)
end

function report(config, grids, rows, assessments, out)
    # Require the candidate AND every tested finer mesh to meet the criteria.
    candidate = findfirst(i -> all(a -> a.eligible, assessments[i:end]), eachindex(assessments))
    empirical = findfirst(i -> all(a -> a.accuracy_pass, assessments[i:end]), eachindex(assessments))
    open(joinpath(out, "study_report.md"), "w") do io
        println(io, "# Normal-shock grid refinement study\n")
        println(io, "Analyzed $(length(grids)) grids. Uniform 1D mesh size is h=L/N (not cbrt(V/N)).\n")
        println(io, "Targets (fractions): field cell-average L1/jump ≤ $(config.field_tol); shock error/L ≤ $(config.shock_tol); scalar GCI/scale ≤ $(config.gci_tol). These are configurable engineering choices, not universal acceptance standards.\n")
        println(io, "| Cells | h (m) | Max normalized L1 | Shock error/L | Accuracy passes | GCI supported |")
        println(io, "|---:|---:|---:|---:|:---:|:---:|")
        for a in assessments
            @printf(io, "| %d | %.6e | %.6e | %.6e | %s | %s |\n", a.cells, a.dx_m,
                    a.max_field_L1_over_jump, a.shock_error_over_length,
                    string(a.accuracy_pass), string(a.gci_supported))
        end
        println(io, "\n## Mesh choice\n")
        if isnothing(candidate)
            println(io, "No tested grid meets all accuracy and supported-GCI criteria across itself and all finer levels. No GCI-qualified mesh recommendation is justified.")
        else
            a = assessments[candidate]
            println(io, "Conditional mesh recommendation: $(a.cells) cells, h=$(a.dx_m) m. It and every tested finer grid meet the configured criteria. Temporal convergence must still be established.")
        end
        if !isnothing(empirical)
            a = assessments[empirical]
            println(io, "\nAccuracy-only candidate: $(a.cells) cells, h=$(a.dx_m) m. This uses known analytical errors and is not a spatial-uncertainty certification.")
        end
        println(io, "\n## Finest-grid scalar uncertainty\n")
        println(io, "All scalar values are signed analytical defects; exact=0. Shock uses domain length; fields use their analytical jump. Physical GCI units are meters for shock and the field's units for mean/plateau defects.\n")
        println(io, "| Quantity | p | GCI/scale | Fourth-grid ratio | Supported | Exact covered |")
        println(io, "|---|---:|---:|---:|:---:|:---:|")
        for r in filter(r -> r.fine_cells==grids[end].n && r.stride==1, rows)
            @printf(io, "| %s | %.4f | %.6e | %.4f | %s | %s |\n", r.quantity,
                    r.observed_order, r.fine_gci_scaled, r.fourth_grid_ratio,
                    string(r.supported), string(r.exact_covered))
        end
        println(io, "\n## Method and limits\n")
        println(io, "Errors are recomputed from profiles, avoiding the old log-parser field-overwrite bug. Both point-sampled and exact cell-average step references are tabulated. Field norms are normalized by the jump, not mixed physical units. Cell averaging of primitive T and u is a diagnostic comparison, not a proof of conservative variable accuracy.")
        println(io, "\nThe exact step location and time come from each verification log. Signed shock and domain-mean defects subtract the moving exact solution at that time. This removes known exact motion but does not establish common-time numerical equivalence. Logs round shock locations to 1e-9 m; smaller errors are not resolved by these inputs.")
        println(io, "\nTriplets use actual unequal spacing ratios. We solve (phi3-phi2)/(phi2-phi1)=(h3^p-h2^p)/(h2^p-h1^p), extrapolate phi0=phi1+(phi1-phi2)/(r21^p-1), and estimate U1=1.25*abs(phi1-phi2)/(r21^p-1). Intervals are phi1 ± U1. GCI is an uncertainty indicator, not a statistical confidence interval or total uncertainty.")
        println(io, "\nOscillating, stagnant, or non-positive-order sequences receive no extrapolation. Local fitted orders above 2 are retained for inspection but excluded from mesh selection. A fitted triplet alone is not evidence of asymptotic convergence: an independent fourth-grid difference must agree with the fitted model within a ratio of 0.8–1.25. This threshold and the order ceiling are conservative study heuristics. Exact-error coverage is reported separately and required for selection.")
        println(io, "\nAdjacent and every-other-grid triplets are reported. The latter increase ratio separation to approximately 1.44 for a 20% sequence. Selection uses adjacent triplets consistently; every-other triplets are sensitivity diagnostics. No pointwise GCI is applied at the discontinuity. Linf and shock-position errors may oscillate with shock/cell alignment; a small value on one grid is insufficient.")
        println(io, "\nAdaptive timesteps vary with mesh. These data contain combined space/time effects; repeat at reduced timesteps or Courant numbers on at least the candidate and finer grids. Compare changes in each quantity with its tolerance and spatial estimate before calling the result a spatial GCI. Model, input, and temporal uncertainty are not quantified or combined here. Mean defects are not conservation-balance residuals; those require boundary flux histories.")
        println(io, "\nFor a follow-up study use systematic grids with ratios around 1.4–2, keep physical inputs fixed, and compare a common saved physical time. Further cell-count predictions from oscillating sequences are deliberately omitted; no reliable required N follows from them.")
        println(io, "\nMethod reference: [NASA NPARC spatial convergence tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).\n")
        println(io, "Outputs: [mesh summary](grid_refinement_summary.csv), [field errors](field_errors.csv), [GCI triplets](gci_triplets.csv), [mesh assessment](mesh_selection.csv).")
    end
    println("Analyzed $(length(grids)) grids; h=$(grids[1].h) to $(grids[end].h) m.")
    println(isnothing(candidate) ? "No GCI-qualified mesh yet." : "Conditional recommended cells: $(assessments[candidate].cells)")
    !isnothing(empirical) && println("Accuracy-only candidate: $(assessments[empirical].cells) cells.")
    println("Report: $(joinpath(out, "study_report.md"))")
end

function draw_plots(result)
    grids, out = result.grids, result.out
    n = [g.n for g in grids]
    panels = []
    for field in FIELDS
        panel = plot(xscale=:log10, yscale=:log10, xlabel="Cells", ylabel="Error / jump",
                     title=string(field), legend=:bottomleft)
        for (k, label) in enumerate(("L1", "L2", "Linf"))
            plot!(panel, n, [g.errors[(field, :normalized_cell_average)][k] for g in grids];
                  label, marker=:circle)
        end
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(2,2), size=(1100,800)), joinpath(out, "grid_refinement_errors.png"))
    thickness = plot(n, [g.thickness/g.h for g in grids]; marker=:circle, label="10–90%",
                     xlabel="Cells", ylabel="Shock thickness (cells)")
    savefig(thickness, joinpath(out, "grid_refinement_thickness.png"))
    shock = plot(n, [g.qois["shock_offset"] for g in grids]; marker=:circle,
                 label="Signed error / L", xlabel="Cells", ylabel="Shock error / domain length")
    savefig(shock, joinpath(out, "grid_refinement_shock.png"))
    panels = []
    for field in FIELDS
        name = "$(field)_plateau_defect"
        rows = filter(r -> r.quantity==name && r.stride==1, result.gci_rows)
        panel = plot(xlabel="Cells", ylabel="Defect / jump", title=string(field), legend=:best)
        plot!(panel, n, [g.qois[name] for g in grids]; label="Post-shock plateau", marker=:circle)
        for supported in (true, false)
            subset = filter(r -> r.supported==supported && isfinite(r.fine_gci_scaled), rows)
            isempty(subset) && continue
            scatter!(panel, [r.fine_cells for r in subset], [r.fine_value_scaled for r in subset];
                     yerror=[r.fine_gci_scaled for r in subset],
                     label=supported ? "Supported GCI" : "Provisional GCI")
        end
        hline!(panel, [0.0]; label="Exact", linestyle=:dash)
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(2,2), size=(1100,800)), joinpath(out, "grid_refinement_gci.png"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = study(options(ARGS))
    using Plots
    draw_plots(result)
end
