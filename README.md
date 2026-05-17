# Hondius Canada Case Analysis

Bayesian classification of the Canadian Andes virus case (onset 14 May 2026) as generation 2 (direct from Hondius case 1) or generation 3 (from the onboard cluster).

## Structure

```
analysis/                          # Canada-specific analysis scripts and report
  canada_case_probability.jl       # Bayesian classification (Julia)
  canada_case_probability_ggplot.R # Figures (R/ggplot2)
  canada_case_secondary_vs_tertiary.toml  # Hypothesis config
  canada_case_generation_report.qmd       # Quarto report
  canada_case_generation_report.bib       # References
andv-linelist-analysis/            # Submodule: epiforecasts/andv-linelist-analysis
  output/posterior.csv             # Epuyén timing posterior (used as input)
external/
  Hondius_hantavirus_h2026/        # Submodule: kraemer-lab/Hondius_hantavirus_h2026
output/                            # Generated outputs (gitignored)
```

## Setup

```bash
git clone --recurse-submodules https://github.com/billyquilty/hondius-canada-case-analysis.git
cd hondius-canada-case-analysis
```

## Running the analysis

**1. Julia classification:**
```bash
julia --project=andv-linelist-analysis analysis/canada_case_probability.jl
```

**2. R figures:**
```r
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
