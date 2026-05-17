# Canada case generation-2 versus generation-3 timing analysis

Case: canada_2026_05_15 (Canada)

Symptom onset: 2026-05-14

Test positive: 2026-05-15

Notes: Working assumptions for an as-yet uncatalogued Canadian case. The generation hypotheses below are evaluated in Julia by combining the fitted transmission-timing posterior for the nominated source cases with the fitted incubation posterior for the Canadian onset date.

## Baseline posterior probabilities

- Generation-2 infection from case 1: 0.088 (95% CrI 0.026 to 0.196); evaluated exposure range 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.912 (95% CrI 0.804 to 0.974); evaluated exposure range 2026-04-17 to 2026-05-13

## Bayesian exposure ranges induced by source-onset timing

- Generation-2 infection from case 1: 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 2026-04-17 to 2026-05-13

## Sensitivity scenarios

### Baseline

- Generation-2 infection from case 1: 0.088 (95% CrI 0.026 to 0.196); prior 1.00; range 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.912 (95% CrI 0.804 to 0.974); prior 1.00; range 2026-04-17 to 2026-05-13

### Source onsets shifted earlier by 1 day

- Generation-2 infection from case 1: 0.072 (95% CrI 0.019 to 0.168); prior 1.00; range 2026-03-31 to 2026-05-12
- Generation-3 infection from generation-2 Hondius cases: 0.928 (95% CrI 0.832 to 0.981); prior 1.00; range 2026-04-16 to 2026-05-12

### Source onsets shifted later by 1 day

- Generation-2 infection from case 1: 0.110 (95% CrI 0.036 to 0.235); prior 1.00; range 2026-04-02 to 2026-05-14
- Generation-3 infection from generation-2 Hondius cases: 0.890 (95% CrI 0.765 to 0.964); prior 1.00; range 2026-04-18 to 2026-05-14

### Generation-2 prior doubled

- Generation-2 infection from case 1: 0.165 (95% CrI 0.052 to 0.334); prior 2.00; range 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.835 (95% CrI 0.666 to 0.948); prior 1.00; range 2026-04-17 to 2026-05-13

### Generation-3 prior doubled

- Generation-2 infection from case 1: 0.046 (95% CrI 0.013 to 0.107); prior 1.00; range 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.954 (95% CrI 0.893 to 0.987); prior 2.00; range 2026-04-17 to 2026-05-13

### Exposure window: 3-day transmission padding

- Generation-2 infection from case 1: 0.086 (95% CrI 0.026 to 0.191); prior 1.00; range 2026-04-03 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.914 (95% CrI 0.809 to 0.974); prior 1.00; range 2026-04-19 to 2026-05-13

### Exposure window: 7-day transmission padding

- Generation-2 infection from case 1: 0.090 (95% CrI 0.027 to 0.200); prior 1.00; range 2026-03-30 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.910 (95% CrI 0.800 to 0.973); prior 1.00; range 2026-04-15 to 2026-05-13

### Generation-3 confirmed sources only

- Generation-2 infection from case 1: 0.084 (95% CrI 0.025 to 0.186); prior 1.00; range 2026-04-01 to 2026-05-13
- Generation-3 infection from generation-2 Hondius cases: 0.916 (95% CrI 0.814 to 0.975); prior 1.00; range 2026-04-17 to 2026-05-13

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
