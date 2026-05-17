# Canada case generation-2 versus generation-3 timing analysis

Case: canada_2026_05_15 (Canada)

Symptom onset: 2026-05-14

Test positive: 2026-05-15

Notes: Working assumptions for an as-yet uncatalogued Canadian case. The generation hypotheses below are evaluated in Julia by combining the fitted transmission-timing posterior for the nominated source cases with the fitted incubation posterior for the Canadian onset date.

## Baseline posterior probabilities

- Generation-2 infection from case 1: 0.261 (95% CrI 0.089 to 0.471); evaluated exposure range 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.739 (95% CrI 0.529 to 0.911); evaluated exposure range 2026-04-17 to 2026-05-10

## Bayesian exposure ranges induced by source-onset timing

- Generation-2 infection from case 1: 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 2026-04-17 to 2026-05-10

## Sensitivity scenarios

### Baseline

- Generation-2 infection from case 1: 0.261 (95% CrI 0.089 to 0.471); prior 1.00; range 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.739 (95% CrI 0.529 to 0.911); prior 1.00; range 2026-04-17 to 2026-05-10

### Source onsets shifted earlier by 1 day

- Generation-2 infection from case 1: 0.219 (95% CrI 0.066 to 0.423); prior 1.00; range 2026-03-31 to 2026-04-10
- Generation-3 infection from generation-2 Hondius cases: 0.781 (95% CrI 0.577 to 0.934); prior 1.00; range 2026-04-16 to 2026-05-09

### Source onsets shifted later by 1 day

- Generation-2 infection from case 1: 0.313 (95% CrI 0.119 to 0.528); prior 1.00; range 2026-04-02 to 2026-04-12
- Generation-3 infection from generation-2 Hondius cases: 0.687 (95% CrI 0.472 to 0.881); prior 1.00; range 2026-04-18 to 2026-05-11

### Generation-2 prior doubled

- Generation-2 infection from case 1: 0.425 (95% CrI 0.169 to 0.654); prior 2.00; range 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.575 (95% CrI 0.346 to 0.831); prior 1.00; range 2026-04-17 to 2026-05-10

### Generation-3 prior doubled

- Generation-2 infection from case 1: 0.146 (95% CrI 0.045 to 0.301); prior 1.00; range 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.854 (95% CrI 0.699 to 0.955); prior 2.00; range 2026-04-17 to 2026-05-10

### Exposure window: 3-day transmission padding

- Generation-2 infection from case 1: 0.285 (95% CrI 0.099 to 0.503); prior 1.00; range 2026-04-03 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.715 (95% CrI 0.497 to 0.901); prior 1.00; range 2026-04-19 to 2026-05-10

### Exposure window: 7-day transmission padding

- Generation-2 infection from case 1: 0.244 (95% CrI 0.082 to 0.449); prior 1.00; range 2026-03-30 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.756 (95% CrI 0.551 to 0.918); prior 1.00; range 2026-04-15 to 2026-05-10

### Generation-3 confirmed sources only

- Generation-2 infection from case 1: 0.250 (95% CrI 0.085 to 0.454); prior 1.00; range 2026-04-01 to 2026-04-11
- Generation-3 infection from generation-2 Hondius cases: 0.750 (95% CrI 0.546 to 0.915); prior 1.00; range 2026-04-17 to 2026-05-10

## Candidate source cases from Hondius line list

- Generation-2 infection from case 1: case 1, status=probable, onset=2026-04-06, confirmation=missing, nationality=dutch
- Generation-3 infection from generation-2 Hondius cases: case 2, status=confirmed, onset=2026-04-24, confirmation=2026-05-04, nationality=dutch
- Generation-3 infection from generation-2 Hondius cases: case 3, status=confirmed, onset=2026-04-22, confirmation=2026-05-02, nationality=british
- Generation-3 infection from generation-2 Hondius cases: case 4, status=confirmed, onset=2026-04-28, confirmation=missing, nationality=german
- Generation-3 infection from generation-2 Hondius cases: case 5, status=confirmed, onset=2026-05-01, confirmation=2026-05-06, nationality=swiss
- Generation-3 infection from generation-2 Hondius cases: case 7, status=confirmed, onset=2026-04-30, confirmation=2026-05-07, nationality=british
- Generation-3 infection from generation-2 Hondius cases: case 8, status=confirmed, onset=2026-04-27, confirmation=2026-05-06, nationality=dutch
- Generation-3 infection from generation-2 Hondius cases: case 12, status=probable, onset=2026-04-28, confirmation=missing, nationality=british
- Generation-3 infection from generation-2 Hondius cases: case 15, status=confirmed, onset=2026-05-10, confirmation=2026-05-11, nationality=french
- Generation-3 infection from generation-2 Hondius cases: case 17, status=monitored, onset=2026-05-10, confirmation=missing, nationality=american
- Generation-3 infection from generation-2 Hondius cases: case 18, status=confirmed, onset=2026-05-10, confirmation=2026-05-11, nationality=spanish

## Sources

- WHO DON600: https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON600
- WHO DON601: https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON601
