# R10 Dorado Backend Validation

Date: 2026-06-09

## Verdict

Icarust R10 simulation is working for the current local Dorado validation path.

The validated claim is:

```text
Icarust can simulate R10 reads from FASTA community inputs.
Dorado can basecall the resulting R10 output.
The basecalled reads map back to the same source references used for simulation.
```

This is a simulation/basecalling/reference-capture claim. It is not a claim that
any adaptive-sampling policy or k-mer target/background decision layer is
production-ready.

## Validated Path

The working path is:

```text
FASTA community manifest
  -> Icarust R10 simulation
  -> FAST5/POD5 conversion path used by the local validation harness
  -> Dorado HAC R10.4.1 E8.2 basecalling
  -> minimap2 mapping back to the simulation source reference
```

The successful local Dorado model used in downstream validation was:

```text
dna_r10.4.1_e8.2_400bps_hac@v5.2.0
```

The compatible sample rate for the current Dorado HAC model bundle is:

```text
5000 Hz
```

Avoid rerunning old 4000 Hz Dorado HAC rows unless a model bundle explicitly
states that it accepts 4000 Hz input.

## Important Code Paths

R10-related files in this repository:

```text
scripts/make_r10_truthset_fast5.py
scripts/analyze_r10_truthset_fastq.py
src/simulation.rs
src/impl_services/data.rs
Profile_tomls/config_dnar10_5khz.toml
Profile_tomls/config_dnar10_5khz_human_barcoded.toml
static/dna_r10.4.1_e8.2_400bps/R10_model.tsv
static/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt
```

Recent implementation commits on `feature/r10-dorado-backend`:

```text
72fdbdc Fix local R10 simulation output path
4ff6f72 Allow custom Squigulator R10 kmer models
85041cd Add Squigulator R10 truth-set backend
359bf80 Add R10 truth-set diagnostics
```

## Validation Evidence

Detailed run reports live in the workflow repository:

```text
/home/shaman.narayanasamy/repositories/adaptive_sampling_microbiome
```

Most relevant documents:

```text
docs/r10_reduced_atcc20_simulation_validation_result.md
docs/r10_reduced_atcc20_generic_manifest_validation_result.md
docs/r10_non_atcc_bacterial_refseq_validation_result.md
docs/r10_final_simulation_generalisability_claim.md
```

The reduced ATCC20 and non-ATCC validation runs demonstrated that simulated R10
reads basecalled by Dorado HAC v5.2.0 map back to the source references used to
produce them.

## Boundaries

Solved:

- R10 simulation can produce usable output for Dorado basecalling.
- Dorado basecalls map back to the source references.
- The generic community-manifest path works for reduced ATCC20 and for a small
  non-ATCC bacterial RefSeq fixture.

Not solved here:

- Calibrated abundance accuracy.
- Live signal/event timing equivalence to a real MinKNOW run.
- Production-scale full ATCC20 validation.
- Adaptive-sampling target/background decision quality.
- k-mer replay transferability across bacterial communities.

The adaptive-sampling replay/controller work belongs in the workflow repository,
not in this Icarust repo.
