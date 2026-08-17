from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from string import Formatter
import tempfile
from typing import Literal

from .fastq import FastqRecord, read_fastq
from .umi import UmiMode, extract_cpm_r2_umi, is_valid_cpm_umi


Orientation = Literal["forward", "reverse-complement"]
ReadSelection = Literal["both", "r1", "r2"]

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def _paths_refer_to_same_location(left: Path, right: Path) -> bool:
    try:
        if left.exists() and right.exists() and left.samefile(right):
            return True
    except OSError:
        pass

    left_key = os.path.normcase(str(left.resolve(strict=False)))
    right_key = os.path.normcase(str(right.resolve(strict=False)))
    return left_key == right_key


def ensure_distinct_paths(**named_paths: str | Path | None) -> None:
    """Reject path aliases before any input is read or output is opened."""

    paths = [(name, Path(value)) for name, value in named_paths.items() if value is not None]
    for index, (left_name, left_path) in enumerate(paths):
        for right_name, right_path in paths[index + 1 :]:
            if _paths_refer_to_same_location(left_path, right_path):
                raise ValueError(
                    f"paths must be distinct: {left_name} and {right_name} "
                    "refer to the same location"
                )


@dataclass(frozen=True)
class ReadTransform:
    orientation: Orientation = "forward"
    trim_left: int = 0
    trim_right: int = 0


@dataclass
class PrepareStats:
    total_pairs: int = 0
    records_written: int = 0
    r1_written: int = 0
    r2_written: int = 0
    skipped_too_short: int = 0
    skipped_n_rate: int = 0
    umi_extracted_pairs: int = 0
    umi_missing_pairs: int = 0


@dataclass(frozen=True)
class _ReadPair:
    read_id: str
    r1: FastqRecord
    r2: FastqRecord
    umi: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def transform_sequence(record: FastqRecord, transform: ReadTransform) -> str:
    if transform.trim_left < 0 or transform.trim_right < 0:
        raise ValueError("trim values must be 0 or greater")

    sequence = record.sequence.upper()
    if transform.trim_left:
        sequence = sequence[transform.trim_left :]
    if transform.trim_right:
        sequence = sequence[: -transform.trim_right]

    if transform.orientation == "forward":
        return sequence
    if transform.orientation == "reverse-complement":
        return reverse_complement(sequence)
    raise ValueError(f"unsupported orientation: {transform.orientation}")


def n_rate(sequence: str) -> float:
    if not sequence:
        return 1.0
    return sequence.count("N") / len(sequence)


def should_write_sequence(sequence: str, *, min_length: int, max_n_rate: float) -> tuple[bool, str | None]:
    if min_length < 0:
        raise ValueError("min_length must be 0 or greater")
    if not 0 <= max_n_rate <= 1:
        raise ValueError("max_n_rate must be between 0 and 1")
    if len(sequence) < min_length:
        return False, "too_short"
    if n_rate(sequence) > max_n_rate:
        return False, "n_rate"
    return True, None


def write_fasta_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def make_query_name(read_id: str, read_label: Literal["R1", "R2"], template: str, umi: str = "") -> str:
    try:
        query_name = template.format(read_id=read_id, read=read_label, umi=umi)
    except (KeyError, IndexError, AttributeError, ValueError) as exc:
        raise ValueError("query name template may only use {read_id}, {read}, and {umi}") from exc
    if not query_name or any(char.isspace() for char in query_name):
        raise ValueError("query names must be non-empty and contain no whitespace")
    return query_name


def _validate_query_name_template_contract(template: str, *, umi_mode: UmiMode) -> None:
    try:
        parsed_template = list(Formatter().parse(template))
    except ValueError as exc:
        raise ValueError(f"invalid query name template: {exc}") from exc

    fields = [field_name for _, field_name, _, _ in parsed_template if field_name is not None]

    allowed_fields = {"read_id", "read", "umi"}
    unsupported_fields = sorted(set(fields) - allowed_fields)
    if unsupported_fields:
        raise ValueError(
            "query name template may only use {read_id}, {read}, and {umi}; "
            f"unsupported field(s): {', '.join(unsupported_fields)}"
        )
    for _, field_name, format_spec, conversion in parsed_template:
        if field_name is not None and (format_spec or conversion is not None):
            raise ValueError(
                "query name template fields may not use format specifications or conversions"
            )
    if fields.count("read_id") != 1:
        if "read_id" not in fields:
            raise ValueError("query name template must include the {read_id} placeholder")
        raise ValueError("query name template must include {read_id} exactly once")
    if fields.count("read") != 1:
        if "read" not in fields:
            raise ValueError("query name template must include the {read} placeholder")
        raise ValueError("query name template must include {read} exactly once")
    if umi_mode == "cpm-r2" and fields.count("umi") != 1:
        if "umi" not in fields:
            raise ValueError("CPM query name template must include the {umi} placeholder")
        raise ValueError("CPM query name template must include {umi} exactly once")
    if fields.count("umi") > 1:
        raise ValueError("query name template may include {umi} at most once")

    template_parts = template.split("|")
    if any(part in {"R1", "R2"} for part in template_parts):
        raise ValueError(
            "query name template may not contain literal pipe-delimited R1/R2 components"
        )
    if "{read}" not in template_parts:
        raise ValueError(
            "query name template must emit {read} as a pipe-delimited R1/R2 component"
        )
    read_component_index = template_parts.index("{read}")
    read_id_component_indices: list[int] = []
    for index, component in enumerate(template_parts):
        try:
            component_fields = [
                field_name
                for _, field_name, _, _ in Formatter().parse(component)
                if field_name is not None
            ]
        except ValueError as exc:
            raise ValueError(f"invalid query name template: {exc}") from exc
        if "read_id" in component_fields:
            read_id_component_indices.append(index)
    if not read_id_component_indices:
        raise ValueError("query name template must include the {read_id} placeholder")
    if max(read_id_component_indices) >= read_component_index:
        raise ValueError(
            "query name template must place the {read} component after {read_id}"
        )
    if umi_mode == "cpm-r2":
        if "UMI={umi}" not in template_parts:
            raise ValueError(
                "CPM query name template must emit UMI={umi} as a pipe-delimited component"
            )
        umi_components = [part for part in template_parts if part.startswith("UMI=")]
        if umi_components != ["UMI={umi}"]:
            raise ValueError(
                "CPM query name template must contain exactly one pipe-delimited "
                "UMI component, and it must be UMI={umi}"
            )
        umi_component_index = template_parts.index("UMI={umi}")
        if umi_component_index <= max(read_id_component_indices):
            raise ValueError(
                "CPM query name template must place UMI={umi} after the {read_id} component"
            )

    probe_umi = "ACGTACGTACGT" if umi_mode == "cpm-r2" else ""
    names = {
        label: make_query_name("PAIR_ID", label, template, probe_umi)
        for label in ("R1", "R2")
    }
    pair_ids: list[str] = []
    for label, name in names.items():
        parts = name.split("|")
        parsed_label: str | None = None
        parsed_index: int | None = None
        for index in range(len(parts) - 1, -1, -1):
            if parts[index] in ("R1", "R2"):
                parsed_label = parts[index]
                parsed_index = index
                break
        if parsed_label != label or parsed_index is None:
            raise ValueError(
                "query name template is parsed ambiguously by the downstream R1/R2 contract"
            )
        pair_ids.append("|".join(parts[:parsed_index] + parts[parsed_index + 1 :]))
    if pair_ids[0] != pair_ids[1]:
        raise ValueError("query name template must produce the same pair ID for R1 and R2")
    if umi_mode == "cpm-r2" and f"UMI={probe_umi}" not in names["R1"].split("|"):
        raise ValueError(
            "CPM query name template must emit UMI={umi} as a pipe-delimited component"
        )


def _write_one_read(
    *,
    fasta,
    stats: PrepareStats,
    read_id: str,
    read_label: Literal["R1", "R2"],
    record: FastqRecord,
    transform: ReadTransform,
    query_name_template: str,
    umi: str,
    min_length: int,
    max_n_rate: float,
) -> None:
    sequence = transform_sequence(record, transform)
    should_write, reason = should_write_sequence(sequence, min_length=min_length, max_n_rate=max_n_rate)
    if not should_write:
        if reason == "too_short":
            stats.skipped_too_short += 1
        elif reason == "n_rate":
            stats.skipped_n_rate += 1
        return

    query_name = make_query_name(read_id, read_label, query_name_template, umi)
    write_fasta_record(fasta, query_name, sequence)
    stats.records_written += 1
    if read_label == "R1":
        stats.r1_written += 1
    else:
        stats.r2_written += 1


def _write_pair(
    *,
    fasta,
    stats: PrepareStats,
    pair: _ReadPair,
    read_selection: ReadSelection,
    r1_transform: ReadTransform,
    r2_transform: ReadTransform,
    query_name_template: str,
    min_length: int,
    max_n_rate: float,
) -> None:
    if read_selection in ("both", "r1"):
        _write_one_read(
            fasta=fasta,
            stats=stats,
            read_id=pair.read_id,
            read_label="R1",
            record=pair.r1,
            transform=r1_transform,
            query_name_template=query_name_template,
            umi=pair.umi,
            min_length=min_length,
            max_n_rate=max_n_rate,
        )
    if read_selection in ("both", "r2"):
        _write_one_read(
            fasta=fasta,
            stats=stats,
            read_id=pair.read_id,
            read_label="R2",
            record=pair.r2,
            transform=r2_transform,
            query_name_template=query_name_template,
            umi=pair.umi,
            min_length=min_length,
            max_n_rate=max_n_rate,
        )


def prepare_paired_fastq_to_fasta(
    r1_path: str | Path,
    r2_path: str | Path,
    fasta_path: str | Path,
    *,
    read_selection: ReadSelection = "both",
    r1_transform: ReadTransform = ReadTransform("forward"),
    r2_transform: ReadTransform = ReadTransform("reverse-complement"),
    query_name_template: str = "{read_id}|{read}",
    min_length: int = 0,
    max_n_rate: float = 1.0,
    strict_ids: bool = True,
    umi_mode: UmiMode = "none",
    umi_anchor_max_mismatches: int = 2,
) -> PrepareStats:
    r1_path = Path(r1_path)
    r2_path = Path(r2_path)
    fasta_path = Path(fasta_path)
    ensure_distinct_paths(r1_path=r1_path, r2_path=r2_path, fasta_path=fasta_path)
    stats = PrepareStats()

    if read_selection not in ("both", "r1", "r2"):
        raise ValueError("read_selection must be one of: both, r1, r2")
    if umi_mode not in ("none", "cpm-r2"):
        raise ValueError("umi_mode must be one of: none, cpm-r2")
    _validate_query_name_template_contract(query_name_template, umi_mode=umi_mode)

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    r1_iter = read_fastq(r1_path, expected_mate=1)
    r2_iter = read_fastq(r2_path, expected_mate=2)

    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{fasta_path.name}.",
        suffix=".tmp",
        dir=fasta_path.parent,
    )
    temporary_path = Path(temporary_name)
    try:
        try:
            temporary_handle = os.fdopen(file_descriptor, "wt", encoding="utf-8", newline="\n")
        except Exception:
            os.close(file_descriptor)
            raise

        with temporary_handle as fasta:
            while True:
                try:
                    r1 = next(r1_iter)
                except StopIteration:
                    try:
                        next(r2_iter)
                    except StopIteration:
                        break
                    raise ValueError("R2 has more records than R1")

                try:
                    r2 = next(r2_iter)
                except StopIteration as exc:
                    raise ValueError("R1 has more records than R2") from exc

                if strict_ids and r1.read_id != r2.read_id:
                    raise ValueError(f"read ID mismatch: R1={r1.read_id!r}, R2={r2.read_id!r}")
                if "|" in r1.read_id or "|" in r2.read_id:
                    raise ValueError("FASTQ read IDs may not contain the reserved '|' delimiter")

                stats.total_pairs += 1
                umi = ""
                if umi_mode == "cpm-r2":
                    extracted_umi = extract_cpm_r2_umi(
                        r2.sequence,
                        max_anchor_mismatches=umi_anchor_max_mismatches,
                    )
                    if is_valid_cpm_umi(extracted_umi):
                        stats.umi_extracted_pairs += 1
                        umi = extracted_umi or "NA"
                    else:
                        stats.umi_missing_pairs += 1
                        # Preserve an ambiguous 12-mer in the query name for
                        # pair-level audit, while classifying it as unusable
                        # for family counting. Structural extraction failure
                        # remains the explicit NA sentinel.
                        umi = extracted_umi or "NA"

                pair = _ReadPair(r1.read_id, r1, r2, umi)
                _write_pair(
                    fasta=fasta,
                    stats=stats,
                    pair=pair,
                    read_selection=read_selection,
                    r1_transform=r1_transform,
                    r2_transform=r2_transform,
                    query_name_template=query_name_template,
                    min_length=min_length,
                    max_n_rate=max_n_rate,
                )

        os.replace(temporary_path, fasta_path)
    finally:
        temporary_path.unlink(missing_ok=True)

    return stats
