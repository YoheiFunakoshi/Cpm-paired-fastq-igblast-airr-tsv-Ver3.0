from __future__ import annotations

from typing import Literal


UmiMode = Literal["none", "cpm-r2"]

CPM_R2_ANCHOR = "TATCAACGCAGAGTGGCCAT"


def is_valid_cpm_umi(umi: str | None) -> bool:
    """Return whether *umi* is an unambiguous CPM 12-mer family key."""

    return bool(umi) and len(umi) == 12 and all(base in "ACGT" for base in umi.upper())


def extract_cpm_r2_umi(
    sequence: str,
    *,
    anchor: str = CPM_R2_ANCHOR,
    max_anchor_mismatches: int = 2,
) -> str | None:
    """Extract the CPM R2 UMI without modifying the query sequence.

    CPM R2 reads are expected to start with:
    anchor + NNNN + T + NNNN + T + NNNN + TCTT + insert.
    The returned UMI is the three NNNN blocks concatenated. Ambiguous bases
    are preserved for audit; :func:`is_valid_cpm_umi` determines whether the
    value is usable as a family key.
    """

    if max_anchor_mismatches < 0:
        raise ValueError("max_anchor_mismatches must be 0 or greater")

    sequence = sequence.upper()
    anchor = anchor.upper()
    umi_end = len(anchor) + 14
    if len(sequence) < umi_end:
        return None

    observed_anchor = sequence[: len(anchor)]
    anchor_mismatches = sum(
        1 for observed, expected in zip(observed_anchor, anchor) if observed != expected
    )
    if anchor_mismatches > max_anchor_mismatches:
        return None

    first = sequence[len(anchor) : len(anchor) + 4]
    second = sequence[len(anchor) + 5 : len(anchor) + 9]
    third = sequence[len(anchor) + 10 : len(anchor) + 14]
    umi = first + second + third
    if len(umi) != 12:
        return None
    return umi
