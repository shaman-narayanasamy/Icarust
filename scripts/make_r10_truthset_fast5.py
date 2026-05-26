#!/usr/bin/env python3
"""Build a tiny R10 FAST5 truth set for Dorado/Icarust debugging.

This script deliberately bypasses the MinKNOW server path and writes controlled
FAST5 reads from Icarust's bundled R10 9-mer model.  The output matrix tests the
main failure axes we care about:

* sequence orientation: forward, reverse, complement, reverse-complement
* raw ADC scaling: ONT inverse calibration versus an intentionally alternate
  formula
* FAST5 metadata calibration: modern R10.4.1 versus legacy R9-like metadata

The resulting FAST5s can be basecalled with Dorado and scored with
``analyze_r10_truthset_fastq.py``.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import uuid
from dataclasses import dataclass
from typing import Iterable

import numpy as np


DNA = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Calibration:
    name: str
    digitisation: float
    offset: float
    range: float
    flow_cell_product_code: str
    exp_script_name: str
    sequencing_kit: str

    @property
    def scale(self) -> float:
        return self.range / self.digitisation


CALIBRATIONS = {
    "modern": Calibration(
        name="modern",
        digitisation=2048.0,
        offset=-243.0,
        range=2048.0 * 0.1462070643901825,
        flow_cell_product_code="FLO-MIN114",
        exp_script_name="sequencing/sequencing_MIN114_DNA,FLO-MIN114,SQK-LSK114",
        sequencing_kit="sqk-lsk114",
    ),
    "legacy": Calibration(
        name="legacy",
        digitisation=8192.0,
        offset=6.0,
        range=1500.0,
        flow_cell_product_code="FLO-MIN106",
        exp_script_name="sequencing/sequencing_MIN106_DNA,FLO-MIN106,SQK-LSK109",
        sequencing_kit="sqk-lsk109",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create controlled R10 FAST5 reads for Dorado diagnostics."
    )
    parser.add_argument(
        "--icarust-root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1],
        help="Icarust repository root containing static/dna_r10.4.1_e8.2_400bps.",
    )
    parser.add_argument(
        "--out-dir",
        type=pathlib.Path,
        required=True,
        help="Directory where FAST5/FASTA/manifest files will be written.",
    )
    parser.add_argument("--sample-rate", type=int, default=5000)
    parser.add_argument("--sequencing-speed", type=int, default=400)
    parser.add_argument("--read-bases", type=int, default=2400)
    parser.add_argument(
        "--noise",
        choices=["none", "laplace"],
        default="none",
        help="Optional deterministic noise. Start with none for diagnosis.",
    )
    parser.add_argument(
        "--model-file",
        type=pathlib.Path,
        default=None,
        help="Override the R10 model TSV. Defaults to Icarust R10_model.tsv.",
    )
    return parser.parse_args()


def build_sequence(length: int) -> str:
    """Create a deterministic non-repetitive sequence with all bases balanced."""
    seed = (
        "ACGTTGCATCCGATGACCTAGGCTAATCGGATCGTACGAT"
        "GCTTACGGAATTCGCCGTAGCTAGTTCGATGGCATCGTTA"
    )
    seq = []
    state = 17
    while len(seq) < length:
        for base in seed:
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF
            idx = (DNA.index(base) + state) % 4
            seq.append(DNA[idx])
            if len(seq) == length:
                break
    return "".join(seq)


def orient(seq: str, mode: str) -> str:
    if mode == "forward":
        return seq
    if mode == "reverse":
        return seq[::-1]
    if mode == "complement":
        return seq.translate(COMPLEMENT)
    if mode == "revcomp":
        return seq.translate(COMPLEMENT)[::-1]
    raise ValueError(mode)


def read_model(path: pathlib.Path) -> dict[str, float]:
    model: dict[str, float] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            model[fields[0].upper()] = float(fields[1])
    if len(model) != 4**9:
        raise ValueError(f"Expected 262144 9-mers in {path}, found {len(model)}")
    return model


def kmers(seq: str, k: int = 9) -> Iterable[str]:
    for i in range(0, len(seq) - k + 1):
        yield seq[i : i + k]


def signal_from_model(
    seq: str,
    model: dict[str, float],
    samples_per_base: int,
    calibration: Calibration,
    scaling: str,
    noise: str,
) -> np.ndarray:
    values = np.array([model[kmer] for kmer in kmers(seq)], dtype=np.float64)
    values = np.repeat(values, samples_per_base)
    if noise == "laplace":
        rng = np.random.default_rng(42)
        values = values + rng.laplace(0.0, 1.0 / math.sqrt(2.0), values.size)

    if scaling == "ont_inverse":
        raw = values / calibration.scale - calibration.offset
    elif scaling == "offset_before_scale":
        raw = (values - calibration.offset) / calibration.scale
    else:
        raise ValueError(scaling)

    return np.rint(raw).clip(np.iinfo(np.int16).min, np.iinfo(np.int16).max).astype(np.int16)


def write_fasta(path: pathlib.Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_fast5(
    path: pathlib.Path,
    read_id: str,
    run_id: str,
    signal: np.ndarray,
    calibration: Calibration,
    sample_rate: int,
) -> None:
    context_tags = {
        "barcoding_enabled": "0",
        "experiment_duration_set": "60",
        "experiment_type": "genomic_dna",
        "local_basecalling": "0",
        "package": "bream4",
        "package_version": "6.3.5",
        "sample_frequency": str(sample_rate),
        "sequencing_kit": calibration.sequencing_kit,
    }
    tracking_id = {
        "asic_id": "817405089",
        "asic_id_eeprom": "5661715",
        "asic_temp": "29.357218",
        "asic_version": "IA02D",
        "auto_update": "0",
        "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/",
        "bream_is_standard": "0",
        "configuration_version": "4.4.13",
        "device_id": "X420",
        "device_type": "gridion",
        "distribution_status": "stable",
        "distribution_version": "21.10.8",
        "exp_script_name": calibration.exp_script_name,
        "exp_script_purpose": "sequencing_run",
        "exp_start_time": "2026-01-01T00:00:00+00:00",
        "flow_cell_id": "TRUTHSET",
        "flow_cell_product_code": calibration.flow_cell_product_code,
        "guppy_version": "5.0.17+99baa5b",
        "heatsink_temp": "34.066406",
        "host_product_code": "GRD-X5B003",
        "host_product_serial_number": "NOTFOUND",
        "hostname": "icarust-truthset",
        "installation_type": "nc",
        "local_firmware_file": "1",
        "operating_system": "linux",
        "protocol_group_id": "r10_truthset",
        "protocol_run_id": "SYNTHETIC_RUN",
        "protocol_start_time": "2026-01-01T00:00:00+00:00",
        "protocols_version": "6.3.5",
        "run_id": run_id,
        "sample_id": "r10_truthset",
        "usb_config": "fx3_1.2.4#fpga_1.2.1#bulk#USB300",
        "version": "4.4.3",
    }
    raw_attrs = {
        "duration": int(signal.size),
        "median_before": 100.0,
        "read_id": read_id,
        "read_number": 1,
        "start_mux": 1,
        "start_time": 1,
        "end_reason": 0,
    }
    channel_info = {
        "digitisation": calibration.digitisation,
        "offset": calibration.offset,
        "range": calibration.range,
        "sampling_rate": float(sample_rate),
        "channel_number": 1,
    }

    try:
        from ont_fast5_api.multi_fast5 import MultiFast5File
    except ModuleNotFoundError:
        write_fast5_h5py(path, read_id, signal, raw_attrs, channel_info, tracking_id, context_tags)
        return

    with MultiFast5File(str(path), "w", driver="core") as multi_f5:
        read = multi_f5.create_empty_read(read_id, run_id)
        read.add_raw_data(signal, attrs=raw_attrs)
        read.add_channel_info(channel_info)
        read.add_tracking_id(tracking_id)
        read.add_context_tags(context_tags)


def write_fast5_h5py(
    path: pathlib.Path,
    read_id: str,
    signal: np.ndarray,
    raw_attrs: dict[str, object],
    channel_info: dict[str, object],
    tracking_id: dict[str, str],
    context_tags: dict[str, str],
) -> None:
    """Write the subset of multi-read FAST5 layout Dorado needs."""
    import h5py

    def add_attrs(group, attrs: dict[str, object]) -> None:
        for key, value in attrs.items():
            group.attrs[key] = value

    with h5py.File(path, "w") as handle:
        read_group = handle.create_group(f"read_{read_id}")
        raw_group = read_group.create_group("Raw")
        raw_group.create_dataset("Signal", data=signal, dtype="<i2", compression="gzip")
        add_attrs(raw_group, raw_attrs)
        add_attrs(read_group.create_group("channel_id"), channel_info)
        add_attrs(read_group.create_group("tracking_id"), tracking_id)
        add_attrs(read_group.create_group("context_tags"), context_tags)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    fast5_dir = out_dir / "fast5"
    npy_dir = out_dir / "npy"
    fast5_dir.mkdir(parents=True, exist_ok=True)
    npy_dir.mkdir(parents=True, exist_ok=True)

    model_path = args.model_file or (
        args.icarust_root / "static/dna_r10.4.1_e8.2_400bps/R10_model.tsv"
    )
    model = read_model(model_path)
    samples_per_base = args.sample_rate // args.sequencing_speed
    if samples_per_base <= 0:
        raise ValueError("sample_rate / sequencing_speed must be at least one sample per base")

    original = build_sequence(args.read_bases)
    records: list[tuple[str, str]] = [("original", original)]
    manifest_rows = []
    run_id = uuid.uuid5(uuid.NAMESPACE_URL, f"icarust-r10-truthset:{out_dir}").hex

    for orientation in ("forward", "reverse", "complement", "revcomp"):
        emitted = orient(original, orientation)
        records.append((orientation, emitted))
        for calibration_name, calibration in CALIBRATIONS.items():
            for scaling in ("ont_inverse", "offset_before_scale"):
                variant = f"{orientation}.{calibration_name}.{scaling}.{args.noise}"
                read_id = str(uuid.uuid5(uuid.NAMESPACE_URL, variant))
                signal = signal_from_model(
                    emitted,
                    model,
                    samples_per_base,
                    calibration,
                    scaling,
                    args.noise,
                )
                npy_path = npy_dir / f"{variant}.npy"
                fast5_path = fast5_dir / f"{variant}.fast5"
                np.save(npy_path, signal)
                write_fast5(fast5_path, read_id, run_id, signal, calibration, args.sample_rate)
                manifest_rows.append(
                    {
                        "variant": variant,
                        "read_id": read_id,
                        "fast5": str(fast5_path),
                        "npy": str(npy_path),
                        "orientation": orientation,
                        "calibration": calibration_name,
                        "scaling": scaling,
                        "noise": args.noise,
                        "sample_rate": args.sample_rate,
                        "sequencing_speed": args.sequencing_speed,
                        "samples_per_base": samples_per_base,
                        "expected_name": orientation,
                        "raw_min": int(signal.min()),
                        "raw_max": int(signal.max()),
                        "raw_mean": float(signal.mean()),
                    }
                )

    write_fasta(out_dir / "expected_sequences.fa", records)
    with (out_dir / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    commands = out_dir / "dorado_commands.sh"
    commands.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "MODEL=${MODEL:-/home/users/snarayanasamy/tools/dorado/models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0}\n"
        "DORADO=${DORADO:-dorado}\n"
        "POD5=${POD5:-pod5}\n"
        f"FAST5_DIR={fast5_dir}\n"
        f"POD5_FILE={out_dir / 'truthset.pod5'}\n"
        f"OUT_FASTQ={out_dir / 'dorado.truthset.fastq'}\n"
        '"$POD5" convert fast5 "$FAST5_DIR" --output "$POD5_FILE" --force-overwrite\n'
        '"$DORADO" basecaller -x "${DORADO_DEVICE:-cpu}" "$MODEL" "$POD5_FILE" '
        '--emit-fastq > "$OUT_FASTQ"\n'
        f"python {pathlib.Path(__file__).with_name('analyze_r10_truthset_fastq.py')} "
        f"--manifest {out_dir / 'manifest.tsv'} "
        f"--expected-fasta {out_dir / 'expected_sequences.fa'} "
        f"--fastq {out_dir / 'dorado.truthset.fastq'} "
        f"--out-tsv {out_dir / 'dorado.truthset.score.tsv'}\n"
    )
    commands.chmod(0o755)

    print(f"Wrote {len(manifest_rows)} FAST5 reads to {fast5_dir}")
    print(f"Manifest: {out_dir / 'manifest.tsv'}")
    print(f"Dorado command helper: {commands}")


if __name__ == "__main__":
    main()
