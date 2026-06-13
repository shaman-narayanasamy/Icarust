# Changelog

## Unreleased

### Added

- Added R10 truth-set diagnostics for tiny FAST5/FASTQ validation runs.
- Added a Squigulator-backed R10 truth-set FAST5 generation path.
- Added support for custom Squigulator R10 k-mer current models in
  `scripts/make_r10_truthset_fast5.py`.
- Added R10 5 kHz example profile TOMLs for Dorado-compatible R10.4.1 E8.2
  validation.
- Added local R10 FAST5/POD5/Dorado validation documentation in
  `docs/r10_dorado_backend_validation.md` so the working simulation claim and
  boundaries are recorded in this repository.

### Fixed

- Fixed local R10 simulation output paths so generated FAST5/POD5/basecall
  artifacts can be found by downstream validation scripts.
- Bumped the MinKNOW API dependency used by readfish-facing code so readfish can
  connect to the simulator in the current local environment.

### Validation

- Confirmed that the repaired R10 simulation path can generate reads that Dorado
  HAC v5.2.0 basecalls and that minimap2 maps back to the source references.
- Confirmed this for a reduced ATCC20 bacterial fixture and for an independent
  non-ATCC bacterial RefSeq fixture.
- The current claim is simulation/basecalling/reference-capture correctness. It
  is not an adaptive-sampling policy-performance claim.
- The successful local Dorado path uses the R10.4.1 E8.2 400 bps HAC v5.2.0
  model family with 5000 Hz simulated input.
- Earlier R10 truth-set backend attempts using Squigulator built-in or custom
  current models did not pass the Dorado tiny truth-set gate; they remain
  provenance for the repair path, not validated production backends.

### Known Follow-Ups

- Add explicit equipment profiles for Flongle-like, MinION-like, and
  PromethION-like simulations. Profiles should own channel count, target yield,
  working pore percentage, duration, sample rate, sequencing speed, flowcell
  metadata, and expected read/base scale. These profiles should be treated as
  simulation approximations, not perfect physical device models.
