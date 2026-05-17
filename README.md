# Hondius Canadian Andes virus case analysis

** Note: analysis up to date as of 17 May 2026 **

Bayesian classification of a Canadian Andes virus case (symptom onset 14 May 2026, confirmed positive 15 May 2026) as generation 2 (direct infection from Hondius case 1) or generation 3 (infection from the onboard secondary cluster).

> **Key result:** Baseline posterior support was **73.9% (95% CrI 52.9–91.1%) for generation 3** versus 26.1% (95% CrI 8.9–47.1%) for generation 2. No sensitivity scenario reversed this ordering. The result indicates additional evidence of human-to-human transmission within the Hondius cluster.

## Approach

**Key assumption:** This analysis assumes that human-to-human (H2H) transmission of Andes virus occurred within the Hondius shipboard cluster — that is, that case 1 infected other passengers (generation 2) via H2H contact, and that generation-2 cases could in turn have infected the Canadian case (generation 3). This assumption is well-supported for Andes virus [@martinez2020; @funk_abbott_2026] but is not independently verified for the Hondius cluster. If the Hondius cluster cases were instead independent zoonotic spillovers, the generation-2/3 distinction would not apply.

Posterior distributions for Andes virus incubation period and transmission timing were taken from Funk and Abbott's Bayesian re-estimation of the 2018–19 Epuyén outbreak ([epiforecasts/andv-linelist-analysis](https://epiforecasts.io/andv-linelist-analysis/dev)) and applied to the Hondius line list. For each posterior draw, the Canadian infection date and candidate source were evaluated jointly: support accumulated where a proposed exposure date was simultaneously consistent with the transmission-timing distribution from the source and with the incubation distribution to Canadian onset (14 May 2026). These contributions were combined with equal priors on generation, then normalised.

Two competing hypotheses were evaluated:

| Hypothesis | Assumed generation | Candidate sources | Exposure window |
|---|---|---|---|
| Generation 2 | Direct from case 1 (onset 6 Apr, died 11 Apr) | 1 | 1 Apr – 11 Apr 2026 |
| Generation 3 | From onboard cluster (cases 2–18) | 10 | 17 Apr – 10 May 2026 |

## Sensitivity analyses

| Scenario | Generation-3 posterior |
|---|---|
| Baseline | 73.9% (52.9–91.1%) |
| Source onsets −1 day | 78.1% (57.7–93.4%) |
| Source onsets +1 day | 68.7% (47.2–88.1%) |
| Generation-2 prior ×2 | 57.5% (34.6–83.1%) |
| Generation-3 prior ×2 | 85.4% (69.9–95.5%) |
| Exposure padding 3 days | 71.5% (49.7–90.1%) |
| Exposure padding 7 days | 75.6% (55.1–91.8%) |
| Confirmed gen-3 sources only | 75.0% (54.6–91.5%) |

## Data and code

- **Hondius line list:** Global.health curated data from the Kraemer Lab — [kraemer-lab/Hondius_hantavirus_h2026](https://github.com/kraemer-lab/Hondius_hantavirus_h2026) (CC BY 4.0), which draws on WHO DON600/601, UKHSA, BBC News, and AP News
- **Epuyén outbreak line list:** Martínez et al. (2020) *NEJM* 383:2230–2241 — [doi:10.1056/NEJMoa2009040](https://doi.org/10.1056/NEJMoa2009040)
- **Timing posterior:** Funk & Abbott (2026) — [epiforecasts.io/andv-linelist-analysis](https://epiforecasts.io/andv-linelist-analysis/dev)
- **Canadian case:** PHAC (16 May 2026) — [canada.ca](https://www.canada.ca/en/public-health/news/2026/05/media-update-on-andes-hantavirus-situation0.html)
- **WHO context:** [DON600](https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON600) · [DON601](https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON601)

## Repository structure

```
analysis/                          # Canada-specific analysis scripts and report
  canada_case_probability.jl       # Bayesian classification (Julia)
  canada_case_probability_ggplot.R # Figures (R/ggplot2)
  canada_case_secondary_vs_tertiary.toml  # Hypothesis config
  canada_case_generation_report.qmd       # Quarto report
  canada_case_generation_report.bib       # References
andv-linelist-analysis/            # Submodule: epiforecasts/andv-linelist-analysis
data/
  posterior.csv                    # Epuyén timing posterior (from Funk & Abbott 2026)
external/
  Hondius_hantavirus_h2026/        # Submodule: kraemer-lab/Hondius_hantavirus_h2026
output/                            # Generated outputs (gitignored)
```

## Setup

```bash
git clone --recurse-submodules https://github.com/bquilty25/hondius-canada-case-analysis.git
cd hondius-canada-case-analysis
```

## Running the analysis

**1. Julia classification:**
```bash
julia --project=andv-linelist-analysis analysis/canada_case_probability.jl
```

**2. R figures:**
```bash
Rscript analysis/canada_case_probability_ggplot.R
```

**3. Render report:**
```bash
quarto render analysis/canada_case_generation_report.qmd
```

## Updating submodules

```bash
# Pull latest upstream andv-linelist-analysis
git submodule update --remote andv-linelist-analysis

# Pull latest Hondius linelist
git submodule update --remote external/Hondius_hantavirus_h2026
```
