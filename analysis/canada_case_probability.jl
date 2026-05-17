using CSV
using DataFrames
using Dates
using Distributions
using Printf
using Statistics
using TOML

const DEFAULT_CONFIG = joinpath(@__DIR__, "canada_case_secondary_vs_tertiary.toml")
const DEFAULT_PADDING_DAYS = 5

function usage()
    println("Usage: julia --project=andv-linelist-analysis analysis/canada_case_probability.jl [config_path] [output_dir]")
    println()
    println("If no arguments are supplied, the default config at analysis/canada_case_secondary_vs_tertiary.toml is used.")
end

function repo_root()
    return normpath(joinpath(@__DIR__, ".."))
end

function resolve_path(root, path)
    return isabspath(path) ? path : normpath(joinpath(root, path))
end

function parse_args(args)
    if any(arg -> arg in ("-h", "--help"), args)
        usage()
        exit(0)
    end
    length(args) <= 2 || error("Expected at most two positional arguments: [config_path] [output_dir]")
    config_path = isempty(args) ? DEFAULT_CONFIG : abspath(args[1])
    output_dir = length(args) >= 2 ? abspath(args[2]) : nothing
    return (; config_path, output_dir)
end

function parse_config(path)
    root = repo_root()
    raw = TOML.parsefile(path)

    case_cfg = raw["case"]
    hypotheses = [
        (; name = string(h["name"]),
            label = string(h["label"]),
            description = string(h["description"]),
            prior_weight = Float64(h["prior_weight"]),
            linked_case_ids = Int.(h["linked_case_ids"]),
            rationale = string(h["rationale"]),
            window_rule = haskey(h, "window_rule") ? h["window_rule"] : nothing)
        for h in raw["hypotheses"]
    ]

    return (
        title = string(raw["title"]),
        posterior_path = resolve_path(root, raw["posterior_path"]),
        hondius_linelist_path = resolve_path(root, raw["hondius_linelist_path"]),
        output_dir = isnothing(get(raw, "output_dir", nothing)) ? nothing : resolve_path(root, raw["output_dir"]),
        case = (; case_id = string(case_cfg["case_id"]),
            country = string(case_cfg["country"]),
            symptom_onset = case_cfg["symptom_onset"],
            test_positive = case_cfg["test_positive"],
            notes = string(case_cfg["notes"])),
        who_events = haskey(raw, "who_events") ? raw["who_events"] : Dict{String, Any}(),
        sources = get(raw, "sources", Any[]),
        hypotheses = hypotheses,
    )
end

parse_optional_date(x) = ismissing(x) || x in ("", "NA") ? missing : Date(x)

function prepare_hondius_linelist(linelist)
    if "symptom_onset" in names(linelist)
        linelist[!, "symptom_onset"] = parse_optional_date.(linelist[!, "symptom_onset"])
    end
    if "confirmation_date" in names(linelist)
        linelist[!, "confirmation_date"] = parse_optional_date.(linelist[!, "confirmation_date"])
    end
    return linelist
end

function resolve_hypothesis_windows(hypotheses, hondius, case_cfg; padding_days = DEFAULT_PADDING_DAYS)
    resolved = NamedTuple[]

    for hypothesis in hypotheses
        anchor_cases = filter(row -> row.Gh_ID in Set(hypothesis.linked_case_ids), hondius)
        nrow(anchor_cases) > 0 || error("No Hondius cases found for hypothesis $(hypothesis.name)")

        dated_anchor_cases = filter(row -> !ismissing(row.symptom_onset), anchor_cases)
        nrow(dated_anchor_cases) > 0 || error("No source onset dates available for hypothesis $(hypothesis.name)")

        source_ids = Int.(dated_anchor_cases[!, "Gh_ID"])
        source_onsets = Date.(dated_anchor_cases[!, "symptom_onset"])
        order = sortperm(source_onsets)
        source_ids = source_ids[order]
        source_onsets = source_onsets[order]
        source_prior_weights = fill(1.0 / length(source_ids), length(source_ids))

        exposure_start = minimum(source_onsets) - Day(padding_days)
        exposure_end = case_cfg.symptom_onset - Day(1)
        exposure_start <= exposure_end || error("Resolved exposure window is invalid for hypothesis $(hypothesis.name)")

        push!(resolved, merge(hypothesis, (; exposure_start, exposure_end, source_ids, source_onsets, source_prior_weights)))
    end

    return resolved
end

function summarise_posterior(posterior)
    required = ["μ_inc", "σ_inc", "μ_δ", "σ_δ"]
    missing_cols = setdiff(required, names(posterior))
    isempty(missing_cols) || error("Posterior CSV is missing required columns: $(join(missing_cols, ", "))")

    keep = select(posterior, required)
    rename!(keep,
        "μ_inc" => :mu_inc,
        "σ_inc" => :sigma_inc,
        "μ_δ" => :mu_delta,
        "σ_δ" => :sigma_delta)
    return keep
end

function onset_probability(onset_date, exposure_date, μ_inc, σ_inc)
    inc_days = Dates.value(onset_date - exposure_date)
    if inc_days < 1
        return 0.0
    end
    dist = LogNormal(μ_inc, σ_inc)
    return cdf(dist, inc_days + 1.0) - cdf(dist, inc_days)
end

function transmission_probability(source_onset, exposure_date, μ_delta, σ_delta)
    delta_days = Dates.value(exposure_date - source_onset)
    dist = Normal(μ_delta, σ_delta)
    return cdf(dist, delta_days + 1.0) - cdf(dist, delta_days)
end

function weighted_quantile(values, weights, probability)
    order = sortperm(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = cumsum(sorted_weights) ./ sum(sorted_weights)
    idx = searchsortedfirst(cumulative, probability)
    return sorted_values[clamp(idx, 1, length(sorted_values))]
end

function weighted_interval_summary(values, weights)
    return (
        mean = sum(values .* weights) / sum(weights),
        lower = weighted_quantile(values, weights, 0.025),
        median = weighted_quantile(values, weights, 0.5),
        upper = weighted_quantile(values, weights, 0.975),
    )
end

function source_date_likelihood_array(hypothesis, case_cfg, posterior)
    exposure_days = collect(hypothesis.exposure_start:Day(1):hypothesis.exposure_end)
    likelihoods = Array{Float64, 3}(undef, nrow(posterior), length(exposure_days), length(hypothesis.source_ids))

    for row in 1:nrow(posterior)
        onset_likelihoods = [
            onset_probability(case_cfg.symptom_onset, exposure_day,
                posterior.mu_inc[row], posterior.sigma_inc[row])
            for exposure_day in exposure_days
        ]

        for (source_idx, source_onset) in enumerate(hypothesis.source_onsets)
            for (date_idx, exposure_day) in enumerate(exposure_days)
                likelihoods[row, date_idx, source_idx] = transmission_probability(
                    source_onset,
                    exposure_day,
                    posterior.mu_delta[row],
                    posterior.sigma_delta[row],
                ) * onset_likelihoods[date_idx]
            end
        end
    end

    return exposure_days, likelihoods
end

function joint_hypothesis_source_date_probabilities(hypotheses, case_cfg, posterior, prior_scale)
    conditional = Dict{String, Tuple{Vector{Date}, Matrix{Float64}}}()
    joint = Dict{String, Tuple{Vector{Date}, Matrix{Float64}}}()
    conditional_source = Dict{String, Tuple{Vector{Int}, Matrix{Float64}}}()
    joint_source = Dict{String, Tuple{Vector{Int}, Matrix{Float64}}}()
    numerators = Dict{String, Matrix{Float64}}()
    source_numerators = Dict{String, Matrix{Float64}}()
    hypothesis_totals = Dict{String, Vector{Float64}}()
    draw_likelihoods = zeros(nrow(posterior))

    for hypothesis in hypotheses
        exposure_days, likelihoods = source_date_likelihood_array(hypothesis, case_cfg, posterior)
        date_priors = fill(1.0 / length(exposure_days), length(exposure_days))
        source_priors = hypothesis.source_prior_weights ./ sum(hypothesis.source_prior_weights)
        date_joint_probs = Matrix{Float64}(undef, nrow(posterior), length(exposure_days))
        date_conditional_probs = Matrix{Float64}(undef, nrow(posterior), length(exposure_days))
        source_joint_probs = Matrix{Float64}(undef, nrow(posterior), length(hypothesis.source_ids))
        source_conditional_probs = Matrix{Float64}(undef, nrow(posterior), length(hypothesis.source_ids))
        hypothesis_weight = prior_scale[hypothesis.name] * hypothesis.prior_weight
        hypothesis_total = zeros(nrow(posterior))

        for row in 1:nrow(posterior)
            for (date_idx, date_prior) in enumerate(date_priors)
                date_joint_probs[row, date_idx] = hypothesis_weight * date_prior * sum(
                    source_priors[source_idx] * likelihoods[row, date_idx, source_idx]
                    for source_idx in eachindex(hypothesis.source_ids)
                )
            end

            for (source_idx, source_prior) in enumerate(source_priors)
                source_joint_probs[row, source_idx] = hypothesis_weight * source_prior * sum(
                    date_priors[date_idx] * likelihoods[row, date_idx, source_idx]
                    for date_idx in eachindex(exposure_days)
                )
            end

            hypothesis_total[row] = sum(date_joint_probs[row, :])
            hypothesis_total[row] > 0 || error("Zero likelihood across exposure dates for hypothesis $(hypothesis.name)")
            date_conditional_probs[row, :] = date_joint_probs[row, :] ./ hypothesis_total[row]
            source_conditional_probs[row, :] = source_joint_probs[row, :] ./ hypothesis_total[row]
        end

        conditional[hypothesis.name] = (exposure_days, date_conditional_probs)
        conditional_source[hypothesis.name] = (hypothesis.source_ids, source_conditional_probs)
        numerators[hypothesis.name] = date_joint_probs
        source_numerators[hypothesis.name] = source_joint_probs
        hypothesis_totals[hypothesis.name] = hypothesis_total
        draw_likelihoods .+= hypothesis_total
    end

    any(==(0.0), draw_likelihoods) && error("At least one posterior draw assigns zero likelihood to all hypothesis-date combinations; revise the exposure windows.")
    draw_weights = draw_likelihoods ./ sum(draw_likelihoods)

    for hypothesis in hypotheses
        exposure_days, _ = conditional[hypothesis.name]
        joint[hypothesis.name] = (exposure_days, numerators[hypothesis.name] ./ reshape(draw_likelihoods, :, 1))
        joint_source[hypothesis.name] = (
            hypothesis.source_ids,
            source_numerators[hypothesis.name] ./ reshape(draw_likelihoods, :, 1),
        )
    end

    return (; joint, conditional, joint_source, conditional_source, draw_weights, hypothesis_totals, draw_likelihoods)
end

function interval_summary(values)
    return (; lower = quantile(values, 0.025), median = quantile(values, 0.5), upper = quantile(values, 0.975))
end

function shift_hypothesis(hypothesis, shift_days)
    return merge(hypothesis, (
        ; exposure_start = hypothesis.exposure_start + Day(shift_days),
        exposure_end = hypothesis.exposure_end + Day(shift_days),
        source_onsets = hypothesis.source_onsets .+ Day(shift_days),
        source_ids = hypothesis.source_ids,
        source_prior_weights = hypothesis.source_prior_weights,
    ))
end

function build_scenarios(hypotheses)
    return [
        (; name = "baseline", label = "Baseline", shift_days = 0,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
        (; name = "window_minus_1d", label = "Source onsets shifted earlier by 1 day", shift_days = -1,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
        (; name = "window_plus_1d", label = "Source onsets shifted later by 1 day", shift_days = 1,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
        (; name = "secondary_prior_x2", label = "Generation-2 prior doubled", shift_days = 0,
            prior_scale = Dict(h.name => (h.name == "secondary" ? 2.0 : 1.0) for h in hypotheses)),
        (; name = "tertiary_prior_x2", label = "Generation-3 prior doubled", shift_days = 0,
            prior_scale = Dict(h.name => (h.name == "tertiary" ? 2.0 : 1.0) for h in hypotheses)),
        (; name = "padding_3d", label = "Exposure window: 3-day transmission padding", shift_days = 0,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
        (; name = "padding_7d", label = "Exposure window: 7-day transmission padding", shift_days = 0,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
        (; name = "confirmed_sources_only", label = "Generation-3 confirmed sources only", shift_days = 0,
            prior_scale = Dict(h.name => 1.0 for h in hypotheses)),
    ]
end

function evaluate_scenario(scenario, hypotheses, case_cfg, posterior; alternative_hypotheses = Dict{String, Any}())
    base_hypotheses = get(alternative_hypotheses, scenario.name, hypotheses)
    shifted = [shift_hypothesis(hypothesis, scenario.shift_days) for hypothesis in base_hypotheses]
    posterior_bundle = joint_hypothesis_source_date_probabilities(shifted, case_cfg, posterior, scenario.prior_scale)
    probabilities = Dict(
        hypothesis.name => posterior_bundle.hypothesis_totals[hypothesis.name] ./ posterior_bundle.draw_likelihoods
        for hypothesis in shifted
    )

    rows = NamedTuple[]
    for hypothesis in shifted
        summary = weighted_interval_summary(probabilities[hypothesis.name], posterior_bundle.draw_weights)
        push!(rows, (
            scenario = scenario.name,
            scenario_label = scenario.label,
            hypothesis = hypothesis.name,
            hypothesis_label = hypothesis.label,
            exposure_start = hypothesis.exposure_start,
            exposure_end = hypothesis.exposure_end,
            prior_weight = scenario.prior_scale[hypothesis.name] * hypothesis.prior_weight,
            probability_mean = summary.mean,
            probability_lower = summary.lower,
            probability_median = summary.median,
            probability_upper = summary.upper,
            mean_weighted_likelihood = sum(posterior_bundle.hypothesis_totals[hypothesis.name]) / length(posterior_bundle.hypothesis_totals[hypothesis.name]),
        ))
    end

    return (; shifted, probabilities, draw_weights = posterior_bundle.draw_weights,
        date_distributions = posterior_bundle.joint,
        conditional_date_distributions = posterior_bundle.conditional,
        source_distributions = posterior_bundle.joint_source,
        conditional_source_distributions = posterior_bundle.conditional_source,
        summary = DataFrame(rows))
end

function summarise_date_distributions(scenario_name, scenario_label, hypotheses, date_distributions, conditional_date_distributions, draw_weights)
    rows = NamedTuple[]
    for hypothesis in hypotheses
        dates, probs = date_distributions[hypothesis.name]
        _, conditional_probs = conditional_date_distributions[hypothesis.name]
        for (idx, exposure_date) in enumerate(dates)
            values = probs[:, idx]
            conditional_values = conditional_probs[:, idx]
            summary = weighted_interval_summary(values, draw_weights)
            conditional_summary = weighted_interval_summary(conditional_values, draw_weights)
            push!(rows, (
                scenario = scenario_name,
                scenario_label = scenario_label,
                hypothesis = hypothesis.name,
                hypothesis_label = hypothesis.label,
                exposure_date = exposure_date,
                probability_mean = summary.mean,
                probability_lower = summary.lower,
                probability_median = summary.median,
                probability_upper = summary.upper,
                conditional_probability_mean = conditional_summary.mean,
                conditional_probability_lower = conditional_summary.lower,
                conditional_probability_median = conditional_summary.median,
                conditional_probability_upper = conditional_summary.upper,
            ))
        end
    end
    return DataFrame(rows)
end

function summarise_source_distributions(scenario_name, scenario_label, hypotheses, source_distributions, conditional_source_distributions, draw_weights)
    rows = NamedTuple[]
    for hypothesis in hypotheses
        source_ids, probs = source_distributions[hypothesis.name]
        _, conditional_probs = conditional_source_distributions[hypothesis.name]
        for (idx, source_id) in enumerate(source_ids)
            values = probs[:, idx]
            conditional_values = conditional_probs[:, idx]
            summary = weighted_interval_summary(values, draw_weights)
            conditional_summary = weighted_interval_summary(conditional_values, draw_weights)
            push!(rows, (
                scenario = scenario_name,
                scenario_label = scenario_label,
                hypothesis = hypothesis.name,
                hypothesis_label = hypothesis.label,
                source_id = source_id,
                probability_mean = summary.mean,
                probability_lower = summary.lower,
                probability_median = summary.median,
                probability_upper = summary.upper,
                conditional_probability_mean = conditional_summary.mean,
                conditional_probability_lower = conditional_summary.lower,
                conditional_probability_median = conditional_summary.median,
                conditional_probability_upper = conditional_summary.upper,
            ))
        end
    end
    return DataFrame(rows)
end

function summarise_resolved_windows(hypotheses)
    return DataFrame((
        hypothesis = h.name,
        hypothesis_label = h.label,
        exposure_start = h.exposure_start,
        exposure_end = h.exposure_end,
        linked_case_ids = join(string.(h.linked_case_ids), ";"),
        rationale = h.rationale,
    ) for h in hypotheses)
end

function supporting_cases(linelist, hypotheses)
    rename_map = Dict(
        Symbol("Gh_ID") => :gh_id,
        :status => :status,
        :symptom_onset => :symptom_onset,
        :confirmation_date => :confirmation_date,
        :nationality => :nationality,
        :linked_information => :linked_information,
        :sources => :sources,
    )

    rows = DataFrame()
    for hypothesis in hypotheses
        ids = Set(hypothesis.linked_case_ids)
        slice = filter(row -> row.Gh_ID in ids, linelist)
        if nrow(slice) == 0
            continue
        end
        keep = select(slice, collect(keys(rename_map)))
        rename!(keep, rename_map)
        keep[!, :hypothesis] .= hypothesis.name
        keep[!, :hypothesis_label] .= hypothesis.label
        append!(rows, keep; cols = :union)
    end
    return rows
end

function write_markdown(path, config, baseline_summary, all_summaries, support_rows)
    open(path, "w") do io
        println(io, "# ", config.title)
        println(io)
        println(io, "Case: ", config.case.case_id, " (", config.case.country, ")")
        println(io)
        println(io, "Symptom onset: ", config.case.symptom_onset)
        println(io)
        println(io, "Test positive: ", config.case.test_positive)
        println(io)
        println(io, "Notes: ", config.case.notes)
        println(io)
        println(io, "## Baseline posterior probabilities")
        println(io)
        for row in eachrow(baseline_summary)
            @printf(io, "- %s: %.3f (95%% CrI %.3f to %.3f); evaluated exposure range %s to %s\n",
                row.hypothesis_label,
                row.probability_median,
                row.probability_lower,
                row.probability_upper,
                row.exposure_start,
                row.exposure_end)
        end
        println(io)
        println(io, "## Bayesian exposure ranges induced by source-onset timing")
        println(io)
        for hypothesis in config.hypotheses
            println(io, "- ", hypothesis.label, ": ", hypothesis.exposure_start, " to ", hypothesis.exposure_end)
        end
        println(io)
        println(io, "## Sensitivity scenarios")
        println(io)
        for scenario in unique(all_summaries.scenario)
            scenario_rows = filter(row -> row.scenario == scenario, all_summaries)
            println(io, "### ", first(scenario_rows.scenario_label))
            println(io)
            for row in eachrow(scenario_rows)
                @printf(io, "- %s: %.3f (95%% CrI %.3f to %.3f); prior %.2f; range %s to %s\n",
                    row.hypothesis_label,
                    row.probability_median,
                    row.probability_lower,
                    row.probability_upper,
                    row.prior_weight,
                    row.exposure_start,
                    row.exposure_end)
            end
            println(io)
        end
        println(io, "## Candidate source cases from Hondius line list")
        println(io)
        for row in eachrow(support_rows)
            println(io, "- ", row.hypothesis_label, ": case ", row.gh_id,
                ", status=", row.status,
                ", onset=", row.symptom_onset,
                ", confirmation=", row.confirmation_date,
                ", nationality=", row.nationality)
        end
        println(io)
        println(io, "## Sources")
        println(io)
        for source in config.sources
            println(io, "- ", source["label"], ": ", source["url"])
        end
    end
end

function main(args)
    cli = parse_args(args)
    config = parse_config(cli.config_path)
    output_dir = isnothing(cli.output_dir) ? something(config.output_dir, joinpath(repo_root(), "output", "canada_case_probability")) : cli.output_dir

    mkpath(output_dir)

    posterior = summarise_posterior(CSV.read(config.posterior_path, DataFrame))
    hondius = prepare_hondius_linelist(CSV.read(config.hondius_linelist_path, DataFrame))
    raw_hypotheses = config.hypotheses
    hypotheses_5d = resolve_hypothesis_windows(raw_hypotheses, hondius, config.case; padding_days = 5)
    hypotheses_3d = resolve_hypothesis_windows(raw_hypotheses, hondius, config.case; padding_days = 3)
    hypotheses_7d = resolve_hypothesis_windows(raw_hypotheses, hondius, config.case; padding_days = 7)
    config = merge(config, (; hypotheses = hypotheses_5d))

    confirmed_gh_ids = Set(
        [row.Gh_ID for row in eachrow(hondius)
         if !ismissing(row.status) && lowercase(strip(string(row.status))) == "confirmed"]
    )
    tertiary_confirmed_ids = let
        tertiary = first(filter(h -> h.name == "tertiary", raw_hypotheses))
        filter(id -> id in confirmed_gh_ids, tertiary.linked_case_ids)
    end
    hypotheses_confirmed = if !isempty(tertiary_confirmed_ids)
        confirmed_raw = map(raw_hypotheses) do h
            h.name == "tertiary" ? merge(h, (; linked_case_ids = tertiary_confirmed_ids)) : h
        end
        try
            resolve_hypothesis_windows(confirmed_raw, hondius, config.case; padding_days = 5)
        catch e
            @warn "Skipping confirmed-only scenario: $e"
            nothing
        end
    else
        @warn "No confirmed gen-3 sources found; skipping confirmed-only scenario"
        nothing
    end

    alternative_hypotheses = Dict{String, Any}(
        "padding_3d" => hypotheses_3d,
        "padding_7d" => hypotheses_7d,
    )
    isnothing(hypotheses_confirmed) || (alternative_hypotheses["confirmed_sources_only"] = hypotheses_confirmed)

    scenarios = build_scenarios(config.hypotheses)
    scenarios = filter(s -> s.name != "confirmed_sources_only" || haskey(alternative_hypotheses, s.name), scenarios)
    scenario_results = [evaluate_scenario(scenario, config.hypotheses, config.case, posterior; alternative_hypotheses) for scenario in scenarios]
    summary = reduce(vcat, [result.summary for result in scenario_results])
    baseline = filter(row -> row.scenario == "baseline", summary)
    support_rows = supporting_cases(hondius, first(scenario_results).shifted)
    date_distribution_summary = summarise_date_distributions(
        first(scenario_results).summary.scenario[1],
        first(scenario_results).summary.scenario_label[1],
        first(scenario_results).shifted,
        first(scenario_results).date_distributions,
        first(scenario_results).conditional_date_distributions,
        first(scenario_results).draw_weights,
    )
    source_distribution_summary = summarise_source_distributions(
        first(scenario_results).summary.scenario[1],
        first(scenario_results).summary.scenario_label[1],
        first(scenario_results).shifted,
        first(scenario_results).source_distributions,
        first(scenario_results).conditional_source_distributions,
        first(scenario_results).draw_weights,
    )
    resolved_windows = summarise_resolved_windows(first(scenario_results).shifted)

    draw_rows = DataFrame(draw = 1:nrow(posterior))
    draw_rows[!, :posterior_draw_weight] = first(scenario_results).draw_weights
    for (hypothesis_name, values) in first(scenario_results).probabilities
        draw_rows[!, Symbol("prob_" * hypothesis_name)] = values
    end
    probability_cols = names(draw_rows, r"^prob_")
    draw_rows[!, :probability_sum] = vec(sum(Matrix(draw_rows[:, probability_cols]), dims = 2))

    CSV.write(joinpath(output_dir, "classification_summary.csv"), summary)
    CSV.write(joinpath(output_dir, "baseline_draws.csv"), draw_rows)
    CSV.write(joinpath(output_dir, "baseline_date_distribution.csv"), date_distribution_summary)
    CSV.write(joinpath(output_dir, "baseline_source_distribution.csv"), source_distribution_summary)
    CSV.write(joinpath(output_dir, "candidate_source_cases.csv"), support_rows)
    CSV.write(joinpath(output_dir, "resolved_hypothesis_windows.csv"), resolved_windows)
    write_markdown(joinpath(output_dir, "classification_summary.md"), config, baseline, summary, support_rows)

    println(config.title)
    for row in eachrow(baseline)
        @printf("%s: %.3f (95%% CrI %.3f to %.3f)\n",
            row.hypothesis_label,
            row.probability_median,
            row.probability_lower,
            row.probability_upper)
    end
    println("Wrote outputs to: ", output_dir)
end

main(ARGS)