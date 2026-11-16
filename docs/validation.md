# Validation

This document lists how the results of pedkit were checked. Every check below is an automated test in `tests/` and runs with `python -m pytest`. The principle throughout is that a number computed by pedkit is compared with the same number obtained by a different route: a hand derivation from first principles, a second exact algorithm that shares no code with the first, unbiased simulation, or direct numerical integration.

## 1. Mendelian inference and Bayes risk

### 1.1. Textbook answers derived by hand

| Case | Hand derivation | Test |
|---|---|---|
| Healthy sib of an AR-affected child, both parents healthy | 2/3 carrier | `test_validation.py::test_cf_healthy_sib_carrier_two_thirds` |
| Child of that sib and a healthy married-in partner, q = 0.02 | 2/3 × 2q/(1+q) × 1/4 = 1/153 | `test_cf_child_risk_exact_fraction` |
| DMD: sister of an affected boy with two healthy sons | prior 1/2, conditional (1/2)² → 1/5 carrier | `test_dmd_one_fifth_and_one_twentieth` |
| Pregnancy of unknown sex of that woman | 1/5 × 1/2 × 1/2 = 1/20 affected | same |
| AD: child of an affected heterozygote, q → 0 | 1/2 affected | `test_ad_child_of_affected_heterozygote_one_half` |
| AD with penetrance f: unaffected child of an affected parent | (1 − f)/(2 − f) carrier | `test_reduced_penetrance_bayes` |
| AR healthy married-in founder | 2q/(1+q) carrier | `test_inference.py::test_exact_classic_carrier_probabilities` |
| Under XR only, sister of an affected boy (test removed) | 1/2 carrier | `test_risk.py::test_mode_averaging` |
| AD descent line of 12 generations, q → 0 | (1/2)¹² carrier | `test_validation.py::test_peeling_has_no_size_limit_descent_line` |

### 1.2. Two exact algorithms

The likelihood P(data | mode), the prior probability of the supporting genotype pool and every individual's posterior genotype marginal are computed by generalised Elston–Stewart peeling (`peeling_likelihood`) and, independently, by enumerating every joint genotype assignment (`enumeration_likelihood`). The two share only the genetic model (`pedkit/inference/modes.py`). On all four Mendelian example pedigrees and all seven modes they agree to a relative error of 10⁻¹⁰ (`test_peeling_equals_enumeration`). The Bayes risk table computed from peeling and from enumeration agrees to 10⁻¹² in risk and mode weights (`test_bayes_exact_equals_enumeration`).

Peeling has no size limit, while enumeration stops at 14 individuals. Two tests check peeling beyond that limit against quantities known in closed form: the 25-individual AD descent line above, and a pedigree in which an unrelated 12-person family is added to the CF example, which must leave the risk of III-1 at exactly 1/153 (`test_unrelated_family_does_not_change_the_risk`).

### 1.3. Monte Carlo against exact

The Monte Carlo likelihood is an unbiased estimate with a reported standard error. For every mode of the mitochondrial example, and for selected modes of the CF example with and without the genotype test and with reduced penetrance and new mutation, the estimate lies within four standard errors of the exact value (`test_monte_carlo_likelihood_within_four_standard_errors`, `test_inference.py`). Importance sampling with a proposal frequency is checked to stay unbiased, to enlarge the supporting pool more than twentyfold and to lower the standard error (`test_importance_sampling_is_unbiased_and_grows_the_pool`). The Monte Carlo Bayes risk of the DMD example is within 0.02 of 1/5 (`test_risk.py::test_monte_carlo_agrees_with_exact`).

### 1.4. Monozygotic twins

MZ twins share one genotype in peeling, enumeration and simulation. Under AD with full penetrance the MZ co-twin of an affected individual is affected with probability 1, while a DZ co-twin has probability 1/2 (`test_mz_twin_of_affected_is_affected_under_full_penetrance`). MZ twins with discordant genotypes give likelihood 0 in both exact methods (`test_mz_twins_discordant_genotypes_are_impossible`).

### 1.5. Validation rule W2

A tested `1/1` genotype in an unaffected individual under AR with full penetrance raises W2; the same data with penetrance 0.9 do not (`test_w2_genotype_contradicts_mode`).

## 2. Multifactorial model

### 2.1. Against bivariate normal integrals

For a pair of relatives with liability correlation r the risk of the second given that the first is affected is P(L₂ > T | L₁ > T), a one-dimensional integral computed by the trapezoidal rule on 200 001 points. The Elston–Stewart risk (hypergeometric polygenic model with Richardson extrapolation) agrees with it to 1 % relative error for:

| Relationship | r | K, h² | Test |
|---|---|---|---|
| Sib of an affected sib, parents of unknown status | h²/2 | 0.01, 0.8 and 0.05, 0.5 | `test_mf.py::test_sib_risk_matches_exact_integral` |
| Child of an affected parent, other parent of unknown status | h²/2 | 0.02, 0.7 | `test_first_degree_relative_risk_bivariate_integral[parent]` |
| Sib of an affected sib, parents of unknown status | h²/2 | 0.02, 0.7 | `test_first_degree_relative_risk_bivariate_integral[sib]` |
| MZ co-twin of an affected twin | h² | 0.01, 0.6 | `test_mz_twin_kinship_and_mf_risk` |

In `test_sib_risk_matches_exact_integral` and in the MZ twin test the MCMC estimate is also within five standard errors (plus 10⁻³) of the integral, and in the former its Gelman–Rubin R̂ is below 1.05.

### 2.2. Against direct multivariate normal simulation

For the eight-person pyloric stenosis example with K = 0.05 and h² = 0.6, liabilities are drawn directly from N(0, Σ) with Σ = h²·2Φ + (1 − h²)·I (ten million draws), the draws that match every observed phenotype are kept, and the share of kept draws with the target above the threshold is counted. Both the Elston–Stewart and the MCMC risks lie within four standard errors of this direct estimate (`test_mf_against_direct_multivariate_normal_simulation`). This test shares no code with either method except the kinship matrix, which is checked separately against the textbook coefficients (`test_mf.py::test_kinship`).

### 2.3. Limits and qualitative properties

With h² = 0 the risk equals the prevalence and the relative risk is 1 (`test_zero_heritability_gives_population_risk`). The risk of an affected individual whose relatives are all unknown equals the prevalence (`test_own_phenotype_is_ignored`). The Carter effect holds: with sex-specific prevalences the son of an affected mother has a higher risk than the son of an affected father, and the relative risk of a daughter exceeds that of a son (`test_carter_effect`).

## 3. File format

Reading and writing a file without changes is byte-identical; unknown columns and metadata are preserved; changing one field rewrites only its line (`test_pedx.py`, section 1). Every ERROR rule E1–E9 is triggered by corrupting one line of the example file and must be reported with the right rule and line number; W3 and W4 are tested the same way. The `.fam` export is checked column by column.

## 4. NSGC drawing conventions

The drawing follows the recommendations of Bennett et al. (2008) and the 2022 revision (Bennett et al. 2022). The table lists the conventions implemented and those not implemented in version 0.1. The SVG output is restricted to the SVG 1.2 Tiny feature set (no clip paths, patterns or markers; `test_svg_is_qt_tiny_compatible`).

| Convention | Status in pedkit 0.1 |
|---|---|
| Square male, circle female, diamond unknown sex | implemented |
| Affected: filled symbol; several conditions: symbol divided into sections | implemented (vertical sections, one per disease; patterns solid, hatched, dotted, quadrant) |
| Unknown phenotype marked with `?` | implemented (pedkit convention inside the section) |
| Carrier not expected to manifest: central dot | implemented for tested `0/1` genotypes only |
| Deceased: diagonal line through the symbol | implemented |
| Age below the symbol, `d.` before the age at death | implemented |
| Proband: arrow with `P`; consultand: arrow without `P` | implemented |
| Pregnancy: `P` inside the symbol | implemented (`planned` and `pregnancy`) |
| Spontaneous abortion: small triangle | implemented (`miscarriage`) |
| Stillbirth: deceased symbol with `SB` | implemented |
| Adoption: brackets around the symbol | implemented (`in`, `out`); the dashed line of descent for adoption is not drawn |
| Relationship line; consanguinity double line; separation two slashes | implemented; consanguinity detected from common ancestors or `##union` |
| Sibship line, twins from one point, MZ connected by a horizontal bar | implemented; `?:` zygosity drawn as DZ |
| Generations in Roman numerals | implemented |
| Legend of conditions | implemented |
| Gender identity label (2022) | implemented as a text label below the symbol; the symbol shape follows the `sex` column, whereas the 2022 revision recommends the symbol to reflect gender identity with sex assigned at birth noted below. Users who follow the 2022 recommendation strictly should enter gender identity in `sex` and assigned sex in `gender`. |
| Termination of pregnancy, ectopic pregnancy, infertility, no children by choice, assisted reproduction and donors, "evaluated" asterisk, multiple individuals of one sex in one symbol, documented evaluation | not implemented |

## 5. What the tests do not cover

The tests check that the computations are correct under the stated models. They do not check that the models fit a particular disease (sections 6.8, 7.6, 8.7 and 12 of the user guide), and they do not include a user study of the viewer. The layout algorithm is checked for valid output and symbol count, not for aesthetic quality on large or highly consanguineous pedigrees.

## References

Bennett, R. L., French, K. S., Resta, R. G. and Doyle, D. L. (2008). Standardized human pedigree nomenclature: update and assessment of the recommendations of the National Society of Genetic Counselors. *Journal of Genetic Counseling* 17, 424–433. https://doi.org/10.1007/s10897-008-9169-9

Bennett, R. L., French, K. S., Resta, R. G. and Austin, J. (2022). Practice resource-focused revision: Standardized pedigree nomenclature update centered on sex and gender inclusivity: A practice resource of the National Society of Genetic Counselors. *Journal of Genetic Counseling* 31, 1238–1248. https://doi.org/10.1002/jgc4.1621
