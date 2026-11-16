# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-15

First public release.

### Added

- PEDX 0.1 format: reader, byte-identical writer, PLINK `.fam` export.
- Validation with nine error rules (E1–E9) and four warning rules (W1–W4),
  each reported with its line number.
- NSGC-standard SVG drawing (SVG 1.2 Tiny), PDF/PNG export and a PyQt5 viewer
  (`pedkit-view`).
- Mode-of-inheritance inference over seven Mendelian modes by exact
  generalised Elston–Stewart peeling (no size limit) and by Monte Carlo gene
  dropping with importance sampling, genotype pools and conflict attribution;
  hierarchical diagnosis tree.
- Bayes risk tables (prior, conditional, joint, posterior) for carrier,
  affected, intermediate, mutant and genotype events, with mode averaging;
  exact, enumeration and Monte Carlo methods.
- Multifactorial liability-threshold risk by Elston–Stewart peeling of the
  hypergeometric polygenic model with Richardson extrapolation and by Gibbs
  sampling MCMC; sex-specific thresholds (Carter effect); five diagnostic
  plots; viewer tab.
- MZ twins share one genotype (Mendelian) and one polygenic value
  (multifactorial).
- Unified command line `pedkit validate|draw|fam|infer|risk|mf`.
- Validation test suite against hand derivations, a second exact algorithm,
  simulation and numerical integration (`docs/validation.md`).
