from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Literal

from .xlsx_writer import write_rows_xlsx
from .umi import is_valid_cpm_umi


ReadLabel = Literal["R1", "R2"]
_REQUIRED_AIRR_COLUMNS = {"sequence_id", "v_call", "j_call", "junction_aa", "productive"}


@dataclass(frozen=True)
class DerivedTsvPaths:
    r1_tsv: Path
    r2_tsv: Path
    integrated_tsv: Path
    counts_tsv: Path
    counts_xlsx: Path
    final_productive_counts_tsv: Path
    final_productive_counts_xlsx: Path
    umi_counts_tsv: Path
    umi_counts_xlsx: Path
    final_productive_umi_counts_tsv: Path
    final_productive_umi_counts_xlsx: Path


@dataclass(frozen=True)
class PairSummaryStats:
    total_pairs: int
    r1_rows: int
    r2_rows: int
    junction_aa_conflicts: int
    included_in_counts: int
    unique_final_clonotypes: int
    final_productive_included_in_counts: int
    unique_final_productive_clonotypes: int


def default_derived_tsv_paths(output_tsv: str | Path) -> DerivedTsvPaths:
    output = Path(output_tsv)
    name = output.name
    lower_name = name.lower()
    if lower_name.endswith(".airr.tsv"):
        sample = name[: -len(".airr.tsv")]
    elif lower_name.endswith(".tsv"):
        sample = name[: -len(".tsv")]
    else:
        sample = output.stem

    return DerivedTsvPaths(
        r1_tsv=output.with_name(f"{sample}.R1.airr.tsv"),
        r2_tsv=output.with_name(f"{sample}.R2.airr.tsv"),
        integrated_tsv=output.with_name(f"{sample}.integrated.tsv"),
        counts_tsv=output.with_name(f"{sample}.integrated_counts.tsv"),
        counts_xlsx=output.with_name(f"{sample}.integrated_counts.xlsx"),
        final_productive_counts_tsv=output.with_name(f"{sample}.final_productive_counts.tsv"),
        final_productive_counts_xlsx=output.with_name(f"{sample}.final_productive_counts.xlsx"),
        umi_counts_tsv=output.with_name(f"{sample}.umi_counts.tsv"),
        umi_counts_xlsx=output.with_name(f"{sample}.umi_counts.xlsx"),
        final_productive_umi_counts_tsv=output.with_name(
            f"{sample}.final_productive_umi_counts.tsv"
        ),
        final_productive_umi_counts_xlsx=output.with_name(
            f"{sample}.final_productive_umi_counts.xlsx"
        ),
    )


def split_and_integrate_airr_tsv(
    input_tsv: str | Path,
    paths: DerivedTsvPaths | None = None,
) -> tuple[DerivedTsvPaths, PairSummaryStats]:
    input_path = Path(input_tsv)
    derived = paths or default_derived_tsv_paths(input_path)
    for path in (
        derived.r1_tsv,
        derived.r2_tsv,
        derived.integrated_tsv,
        derived.counts_tsv,
        derived.counts_xlsx,
        derived.final_productive_counts_tsv,
        derived.final_productive_counts_xlsx,
        derived.umi_counts_tsv,
        derived.umi_counts_xlsx,
        derived.final_productive_umi_counts_tsv,
        derived.final_productive_umi_counts_xlsx,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)

    pairs: dict[str, dict[ReadLabel, dict[str, str]]] = {}
    r1_rows = 0
    r2_rows = 0

    with input_path.open("rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            _write_empty_tsv(derived.r1_tsv)
            _write_empty_tsv(derived.r2_tsv)
            _write_integrated_tsv(derived.integrated_tsv, [])
            _write_counts_tsv(derived.counts_tsv, [])
            write_rows_xlsx(
                derived.counts_xlsx,
                COUNTS_FIELDNAMES,
                [],
                sheet_name="integrated_counts",
                header_labels=COUNTS_XLSX_HEADER_LABELS,
            )
            _write_counts_tsv(derived.final_productive_counts_tsv, [])
            write_rows_xlsx(
                derived.final_productive_counts_xlsx,
                COUNTS_FIELDNAMES,
                [],
                sheet_name="final_productive_counts",
                header_labels=COUNTS_XLSX_HEADER_LABELS,
            )
            _write_umi_counts_tsv(derived.umi_counts_tsv, [])
            write_rows_xlsx(
                derived.umi_counts_xlsx,
                UMI_COUNTS_FIELDNAMES,
                [],
                sheet_name="umi_counts",
                header_labels=COUNTS_XLSX_HEADER_LABELS,
                percentage_fields=UMI_COUNTS_XLSX_PERCENTAGE_FIELDS,
            )
            _write_umi_counts_tsv(derived.final_productive_umi_counts_tsv, [])
            write_rows_xlsx(
                derived.final_productive_umi_counts_xlsx,
                UMI_COUNTS_FIELDNAMES,
                [],
                sheet_name="final_productive_umi_counts",
                header_labels=COUNTS_XLSX_HEADER_LABELS,
                percentage_fields=UMI_COUNTS_XLSX_PERCENTAGE_FIELDS,
            )
            return derived, PairSummaryStats(0, 0, 0, 0, 0, 0, 0, 0)

        missing_columns = sorted(_REQUIRED_AIRR_COLUMNS - set(reader.fieldnames))
        if missing_columns:
            raise ValueError(
                "AIRR TSV is missing required column(s): " + ", ".join(missing_columns)
            )

        with (
            derived.r1_tsv.open("wt", encoding="utf-8", newline="") as r1_handle,
            derived.r2_tsv.open("wt", encoding="utf-8", newline="") as r2_handle,
        ):
            r1_writer = csv.DictWriter(r1_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            r2_writer = csv.DictWriter(r2_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            r1_writer.writeheader()
            r2_writer.writeheader()

            for row in reader:
                pair_id, read_label = pair_id_and_read_label(row.get("sequence_id", ""))
                if not pair_id or read_label is None:
                    raise ValueError(
                        "AIRR sequence_id does not satisfy the R1/R2 query-name contract: "
                        f"{row.get('sequence_id', '')!r}"
                    )
                pair_rows = pairs.setdefault(pair_id, {})
                if read_label in pair_rows:
                    raise ValueError(
                        f"AIRR TSV contains duplicate {read_label} row for pair ID {pair_id!r}"
                    )
                pair_rows[read_label] = row
                if read_label == "R1":
                    r1_writer.writerow(row)
                    r1_rows += 1
                else:
                    r2_writer.writerow(row)
                    r2_rows += 1

    integrated_rows = [_integrated_row(pair_id, pair_rows) for pair_id, pair_rows in sorted(pairs.items())]
    _write_integrated_tsv(derived.integrated_tsv, integrated_rows)
    counts_rows = _counts_rows(integrated_rows)
    _write_counts_tsv(derived.counts_tsv, counts_rows)
    write_rows_xlsx(
        derived.counts_xlsx,
        COUNTS_FIELDNAMES,
        counts_rows,
        sheet_name="integrated_counts",
        header_labels=COUNTS_XLSX_HEADER_LABELS,
    )
    final_productive_counts_rows = _counts_rows(integrated_rows, final_productive_only=True)
    _write_counts_tsv(derived.final_productive_counts_tsv, final_productive_counts_rows)
    write_rows_xlsx(
        derived.final_productive_counts_xlsx,
        COUNTS_FIELDNAMES,
        final_productive_counts_rows,
        sheet_name="final_productive_counts",
        header_labels=COUNTS_XLSX_HEADER_LABELS,
    )
    umi_counts_rows = _umi_counts_rows(integrated_rows)
    _write_umi_counts_tsv(derived.umi_counts_tsv, umi_counts_rows)
    write_rows_xlsx(
        derived.umi_counts_xlsx,
        UMI_COUNTS_FIELDNAMES,
        umi_counts_rows,
        sheet_name="umi_counts",
        header_labels=COUNTS_XLSX_HEADER_LABELS,
        percentage_fields=UMI_COUNTS_XLSX_PERCENTAGE_FIELDS,
    )
    final_productive_umi_counts_rows = _umi_counts_rows(
        integrated_rows,
        final_productive_only=True,
    )
    _write_umi_counts_tsv(
        derived.final_productive_umi_counts_tsv,
        final_productive_umi_counts_rows,
    )
    write_rows_xlsx(
        derived.final_productive_umi_counts_xlsx,
        UMI_COUNTS_FIELDNAMES,
        final_productive_umi_counts_rows,
        sheet_name="final_productive_umi_counts",
        header_labels=COUNTS_XLSX_HEADER_LABELS,
        percentage_fields=UMI_COUNTS_XLSX_PERCENTAGE_FIELDS,
    )
    conflicts = sum(1 for row in integrated_rows if row["junction_aa_status"] == "conflict")
    included = sum(1 for row in integrated_rows if row["include_in_counts"] == "true")
    final_productive_included = sum(
        1 for row in integrated_rows if row["include_in_counts"] == "true" and _is_productive(row["final_productive"])
    )
    return derived, PairSummaryStats(
        len(integrated_rows),
        r1_rows,
        r2_rows,
        conflicts,
        included,
        len(counts_rows),
        final_productive_included,
        len(final_productive_counts_rows),
    )


def pair_id_and_read_label(sequence_id: str) -> tuple[str, ReadLabel | None]:
    text = sequence_id.strip()
    if not text or "|" not in text:
        return text, None
    parts = text.split("|")
    for index in range(len(parts) - 1, -1, -1):
        if parts[index] in ("R1", "R2"):
            pair_id = "|".join(parts[:index] + parts[index + 1 :])
            return pair_id, parts[index]  # type: ignore[return-value]
    return text, None


def _write_empty_tsv(path: Path) -> None:
    path.write_text("", encoding="utf-8", newline="")


def _is_value(value: str | None) -> bool:
    if value is None:
        return False
    text = value.strip()
    return bool(text) and text.lower() not in {"na", "n/a", "none", "null"}


_ALLELE_SUFFIX = re.compile(r"\*\d+(?:[A-Z])?$")
_UMI_PATTERN = re.compile(r"(?:^|\|)UMI=([^|\s]+)")


def gene_candidate_set(call: str) -> str:
    genes = set()
    for part in (call or "").split(","):
        gene = _ALLELE_SUFFIX.sub("", part.strip())
        if gene and gene.lower() not in {"na", "n/a", "none", "null", "x"}:
            genes.add(gene)
    return ",".join(sorted(genes))


def _is_productive(value: str) -> bool:
    return value.strip().lower() in {"t", "true", "yes", "1"}


def _format_percent(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.000000%"
    return f"{(numerator / denominator) * 100:.6f}%"


def _junction_aa_exclude_reasons(junction_aa: str) -> list[str]:
    if not _is_value(junction_aa):
        return ["missing_junction_aa"]

    reasons = []
    if "*" in junction_aa:
        reasons.append("junction_aa_has_stop")
    if not junction_aa.startswith("C"):
        reasons.append("junction_aa_not_c_start")
    if not junction_aa.endswith(("W", "F")):
        reasons.append("junction_aa_not_wf_end")
    if not 5 <= len(junction_aa) <= 40:
        reasons.append("junction_aa_length_out_of_range")
    return reasons


def _is_canonical_junction_aa(junction_aa: str) -> bool:
    return not _junction_aa_exclude_reasons(junction_aa)


def _count_inclusion(
    *,
    final_v_call: str,
    final_j_call: str,
    final_junction_aa: str,
) -> tuple[str, str, str, str]:
    unique_v_gene_set = gene_candidate_set(final_v_call)
    unique_j_gene_set = gene_candidate_set(final_j_call)
    reasons = []
    if not unique_v_gene_set:
        reasons.append("missing_v_call")
    if not unique_j_gene_set:
        reasons.append("missing_j_call")
    reasons.extend(_junction_aa_exclude_reasons(final_junction_aa))
    include = not reasons
    return (
        "true" if include else "false",
        "" if include else ";".join(reasons),
        unique_v_gene_set,
        unique_j_gene_set,
    )


def _get(row: dict[str, str] | None, field: str) -> str:
    if row is None:
        return ""
    return row.get(field, "") or ""


def extract_umi(*values: str) -> str:
    for value in values:
        matches = list(_UMI_PATTERN.finditer(value or ""))
        if not matches:
            continue
        # The generated query-name contract puts the authoritative UMI after
        # the read ID. A read ID itself may legitimately begin with ``UMI=``;
        # only the last pipe-delimited UMI component in each value is therefore
        # interpreted as the CPM UMI.
        umi = matches[-1].group(1).strip()
        if _is_value(umi):
            return umi
    return ""


def _is_known_umi(value: str) -> bool:
    return is_valid_cpm_umi(value.strip())


def _choose_junction_aa(r1: str, r2: str) -> tuple[str, str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has:
        if r1 == r2:
            return r1, "match", "both", "same_junction_aa"
        r1_canonical = _is_canonical_junction_aa(r1)
        r2_canonical = _is_canonical_junction_aa(r2)
        if r1_canonical and not r2_canonical:
            return r1, "conflict", "R1", "conflict_canonical_r1"
        if r2_canonical and not r1_canonical:
            return r2, "conflict", "R2", "conflict_canonical_r2"
        return r1, "conflict", "R1", "conflict_r1_priority"
    if r1_has:
        return r1, "r1_only", "R1", "only_r1_has_junction_aa"
    if r2_has:
        return r2, "r2_only", "R2", "only_r2_has_junction_aa"
    return "", "none", "", "no_junction_aa"


def _choose_prefer_r2(r1: str, r2: str) -> tuple[str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has and r1 == r2:
        return r1, "both", "match"
    if r2_has:
        return r2, "R2", "r2_priority"
    if r1_has:
        return r1, "R1", "r1_only"
    return "", "", "none"


def _choose_by_preferred_read(
    r1: str,
    r2: str,
    preferred_read: str,
    *,
    fallback_r2: bool = False,
) -> tuple[str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has and r1 == r2:
        return r1, "both", "match"
    if preferred_read == "both":
        if r1_has:
            return r1, "R1", "both_r1_priority"
        if r2_has:
            return r2, "R2", "both_fallback_r2"
    if preferred_read == "R1" and r1_has:
        return r1, "R1", "preferred_read"
    if preferred_read == "R2" and r2_has:
        return r2, "R2", "preferred_read"
    if fallback_r2 and r2_has:
        return r2, "R2", "fallback_r2"
    if r1_has:
        return r1, "R1", "fallback_r1"
    if r2_has:
        return r2, "R2", "fallback_r2"
    return "", "", "none"


def _integrated_row(pair_id: str, pair_rows: dict[ReadLabel, dict[str, str]]) -> dict[str, str]:
    r1 = pair_rows.get("R1")
    r2 = pair_rows.get("R2")
    r1_junction_aa = _get(r1, "junction_aa")
    r2_junction_aa = _get(r2, "junction_aa")
    final_junction_aa, junction_status, preferred_read, decision_reason = _choose_junction_aa(
        r1_junction_aa,
        r2_junction_aa,
    )

    final_v, v_source, v_reason = _choose_prefer_r2(_get(r1, "v_call"), _get(r2, "v_call"))
    final_j, j_source, j_reason = _choose_by_preferred_read(
        _get(r1, "j_call"),
        _get(r2, "j_call"),
        preferred_read,
        fallback_r2=True,
    )
    final_d, d_source, d_reason = _choose_by_preferred_read(
        _get(r1, "d_call"),
        _get(r2, "d_call"),
        preferred_read,
        fallback_r2=True,
    )
    final_productive, productive_source, productive_reason = _choose_by_preferred_read(
        _get(r1, "productive"),
        _get(r2, "productive"),
        preferred_read,
        fallback_r2=True,
    )
    final_junction, junction_source, junction_reason = _choose_by_preferred_read(
        _get(r1, "junction"),
        _get(r2, "junction"),
        preferred_read,
        fallback_r2=True,
    )

    include_in_counts, exclude_reason, unique_v_gene_set, unique_j_gene_set = _count_inclusion(
        final_v_call=final_v,
        final_j_call=final_j,
        final_junction_aa=final_junction_aa,
    )

    return {
        "pair_id": pair_id,
        "umi": extract_umi(_get(r2, "sequence_id"), _get(r1, "sequence_id"), pair_id),
        "r1_sequence_id": _get(r1, "sequence_id"),
        "r2_sequence_id": _get(r2, "sequence_id"),
        "final_junction_aa": final_junction_aa,
        "junction_aa_status": junction_status,
        "preferred_read": preferred_read,
        "junction_aa_decision_reason": decision_reason,
        "r1_junction_aa": r1_junction_aa,
        "r2_junction_aa": r2_junction_aa,
        "final_junction": final_junction,
        "junction_source": junction_source,
        "junction_decision_reason": junction_reason,
        "r1_junction": _get(r1, "junction"),
        "r2_junction": _get(r2, "junction"),
        "final_v_call": final_v,
        "v_call_source": v_source,
        "v_call_decision_reason": v_reason,
        "r1_v_call": _get(r1, "v_call"),
        "r2_v_call": _get(r2, "v_call"),
        "final_d_call": final_d,
        "d_call_source": d_source,
        "d_call_decision_reason": d_reason,
        "r1_d_call": _get(r1, "d_call"),
        "r2_d_call": _get(r2, "d_call"),
        "final_j_call": final_j,
        "j_call_source": j_source,
        "j_call_decision_reason": j_reason,
        "r1_j_call": _get(r1, "j_call"),
        "r2_j_call": _get(r2, "j_call"),
        "final_productive": final_productive,
        "productive_source": productive_source,
        "productive_decision_reason": productive_reason,
        "r1_productive": _get(r1, "productive"),
        "r2_productive": _get(r2, "productive"),
        "unique_v_gene_set": unique_v_gene_set,
        "unique_j_gene_set": unique_j_gene_set,
        "include_in_counts": include_in_counts,
        "exclude_reason": exclude_reason,
    }


INTEGRATED_FIELDNAMES = [
    "pair_id",
    "umi",
    "r1_sequence_id",
    "r2_sequence_id",
    "final_junction_aa",
    "junction_aa_status",
    "preferred_read",
    "junction_aa_decision_reason",
    "r1_junction_aa",
    "r2_junction_aa",
    "final_junction",
    "junction_source",
    "junction_decision_reason",
    "r1_junction",
    "r2_junction",
    "final_v_call",
    "v_call_source",
    "v_call_decision_reason",
    "r1_v_call",
    "r2_v_call",
    "final_d_call",
    "d_call_source",
    "d_call_decision_reason",
    "r1_d_call",
    "r2_d_call",
    "final_j_call",
    "j_call_source",
    "j_call_decision_reason",
    "r1_j_call",
    "r2_j_call",
    "final_productive",
    "productive_source",
    "productive_decision_reason",
    "r1_productive",
    "r2_productive",
    "unique_v_gene_set",
    "unique_j_gene_set",
    "include_in_counts",
    "exclude_reason",
]


def _write_integrated_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INTEGRATED_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


COUNTS_FIELDNAMES = [
    "unique_v_gene_set",
    "unique_j_gene_set",
    "final_junction_aa",
    "read_pair_count",
    "match_count",
    "conflict_count",
    "r1_only_count",
    "r2_only_count",
    "productive_true_count",
    "canonical_junction_aa_count",
]


UMI_COUNTS_FIELDNAMES = [
    "unique_v_gene_set",
    "unique_j_gene_set",
    "final_junction_aa",
    "umi_family_count",
    "read_pair_count",
    "umi_known_read_pair_count",
    "umi_missing_read_pair_count",
    "inclusive_support_count",
    "inclusive_support_percent",
    "match_count",
    "conflict_count",
    "r1_only_count",
    "r2_only_count",
    "productive_true_count",
    "canonical_junction_aa_count",
]


COUNTS_XLSX_HEADER_LABELS = {
    "final_junction_aa": "final_junction_aa (canonical)",
}

UMI_COUNTS_XLSX_PERCENTAGE_FIELDS = {"inclusive_support_percent"}


def _counts_rows(
    integrated_rows: list[dict[str, str]],
    *,
    final_productive_only: bool = False,
) -> list[dict[str, str]]:
    counts: dict[tuple[str, str, str], dict[str, int | str]] = {}
    for row in integrated_rows:
        if row.get("include_in_counts", "").strip().lower() != "true":
            continue
        if final_productive_only and not _is_productive(row.get("final_productive", "")):
            continue
        key = (
            row.get("unique_v_gene_set", ""),
            row.get("unique_j_gene_set", ""),
            row.get("final_junction_aa", ""),
        )
        bucket = counts.setdefault(
            key,
            {
                "unique_v_gene_set": key[0],
                "unique_j_gene_set": key[1],
                "final_junction_aa": key[2],
                "read_pair_count": 0,
                "match_count": 0,
                "conflict_count": 0,
                "r1_only_count": 0,
                "r2_only_count": 0,
                "productive_true_count": 0,
                "canonical_junction_aa_count": 0,
            },
        )
        bucket["read_pair_count"] = int(bucket["read_pair_count"]) + 1
        status_key = row.get("junction_aa_status", "none")
        if status_key not in {"match", "conflict", "r1_only", "r2_only"}:
            continue
        bucket[f"{status_key}_count"] = int(bucket[f"{status_key}_count"]) + 1

        productive = row.get("final_productive", "").strip().lower()
        if productive in {"t", "true", "yes", "1"}:
            bucket["productive_true_count"] = int(bucket["productive_true_count"]) + 1
        bucket["canonical_junction_aa_count"] = int(bucket["canonical_junction_aa_count"]) + 1

    sorted_rows = sorted(
        counts.values(),
        key=lambda item: (
            -int(item["read_pair_count"]),
            str(item["unique_v_gene_set"]),
            str(item["unique_j_gene_set"]),
            str(item["final_junction_aa"]),
        ),
    )
    return [{field: str(row[field]) for field in COUNTS_FIELDNAMES} for row in sorted_rows]


def _umi_counts_rows(
    integrated_rows: list[dict[str, str]],
    *,
    final_productive_only: bool = False,
) -> list[dict[str, str]]:
    """Count exact raw 12-mer UMIs independently within each BCR key."""

    counts: dict[tuple[str, str, str], dict[str, int | str | set[str]]] = {}
    for row in integrated_rows:
        if row.get("include_in_counts", "").strip().lower() != "true":
            continue
        if final_productive_only and not _is_productive(row.get("final_productive", "")):
            continue

        key = (
            row.get("unique_v_gene_set", ""),
            row.get("unique_j_gene_set", ""),
            row.get("final_junction_aa", ""),
        )
        bucket = counts.setdefault(
            key,
            {
                "unique_v_gene_set": key[0],
                "unique_j_gene_set": key[1],
                "final_junction_aa": key[2],
                "umi_family_count": 0,
                "read_pair_count": 0,
                "umi_known_read_pair_count": 0,
                "umi_missing_read_pair_count": 0,
                "inclusive_support_count": 0,
                "inclusive_support_percent": "0.000000%",
                "_umi_families": set(),
                "match_count": 0,
                "conflict_count": 0,
                "r1_only_count": 0,
                "r2_only_count": 0,
                "productive_true_count": 0,
                "canonical_junction_aa_count": 0,
            },
        )
        bucket["read_pair_count"] = int(bucket["read_pair_count"]) + 1

        umi = row.get("umi", "")
        if _is_known_umi(umi):
            bucket["umi_known_read_pair_count"] = int(bucket["umi_known_read_pair_count"]) + 1
            umi_families = bucket["_umi_families"]
            if isinstance(umi_families, set):
                umi_families.add(umi)
        else:
            bucket["umi_missing_read_pair_count"] = int(bucket["umi_missing_read_pair_count"]) + 1

        status_key = row.get("junction_aa_status", "none")
        if status_key in {"match", "conflict", "r1_only", "r2_only"}:
            count_field = f"{status_key}_count"
            bucket[count_field] = int(bucket[count_field]) + 1

        if _is_productive(row.get("final_productive", "")):
            bucket["productive_true_count"] = int(bucket["productive_true_count"]) + 1
        bucket["canonical_junction_aa_count"] = int(bucket["canonical_junction_aa_count"]) + 1

    for bucket in counts.values():
        umi_families = bucket["_umi_families"]
        bucket["umi_family_count"] = len(umi_families) if isinstance(umi_families, set) else 0
        bucket["inclusive_support_count"] = int(bucket["umi_family_count"]) + int(
            bucket["umi_missing_read_pair_count"]
        )

    sorted_rows = sorted(
        counts.values(),
        key=lambda item: (
            -int(item["inclusive_support_count"]),
            -int(item["umi_family_count"]),
            -int(item["read_pair_count"]),
            str(item["unique_v_gene_set"]),
            str(item["unique_j_gene_set"]),
            str(item["final_junction_aa"]),
        ),
    )
    inclusive_total = sum(int(row["inclusive_support_count"]) for row in sorted_rows)
    output_rows = []
    for row in sorted_rows:
        row["inclusive_support_percent"] = _format_percent(
            int(row["inclusive_support_count"]),
            inclusive_total,
        )
        output_rows.append({field: str(row[field]) for field in UMI_COUNTS_FIELDNAMES})
    return output_rows


def _write_counts_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COUNTS_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_umi_counts_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=UMI_COUNTS_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
