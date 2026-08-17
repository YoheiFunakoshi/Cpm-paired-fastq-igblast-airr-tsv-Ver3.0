from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal, TextIO
import gzip


MateNumber = Literal[1, 2]


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    header: str
    sequence: str
    quality: str


def normalize_read_id(header: str) -> str:
    text = header.strip()
    if text.startswith("@"):
        text = text[1:]
    read_id = text.split()[0]
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]
    return read_id


def mate_number_from_header(header: str) -> MateNumber | None:
    """Return an explicit R1/R2 marker from a FASTQ header when present."""

    text = header.strip()
    if text.startswith("@"):
        text = text[1:]
    parts = text.split()
    if not parts:
        return None

    observed: list[MateNumber] = []
    if parts[0].endswith("/1"):
        observed.append(1)
    elif parts[0].endswith("/2"):
        observed.append(2)

    if len(parts) > 1 and len(parts[1]) > 1 and parts[1][1] == ":":
        if parts[1][0] == "1":
            observed.append(1)
        elif parts[1][0] == "2":
            observed.append(2)

    if not observed:
        return None
    if any(value != observed[0] for value in observed[1:]):
        raise ValueError("FASTQ header has conflicting mate markers")
    return observed[0]


def open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def read_fastq(
    path: str | Path,
    *,
    expected_mate: MateNumber | None = None,
) -> Iterator[FastqRecord]:
    path = Path(path)
    if expected_mate not in (None, 1, 2):
        raise ValueError("expected_mate must be 1, 2, or None")

    with open_text(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                return

            record_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(f"{path}: incomplete FASTQ record at record {record_number}")

            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            if not header.startswith("@"):
                raise ValueError(f"{path}: FASTQ header does not start with @ at record {record_number}")
            if not plus.startswith("+"):
                raise ValueError(f"{path}: FASTQ plus line does not start with + at record {record_number}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: sequence and quality lengths differ at record {record_number}"
                )

            try:
                observed_mate = mate_number_from_header(header)
            except ValueError as exc:
                raise ValueError(f"{path}: {exc} at record {record_number}") from exc
            if expected_mate is not None and observed_mate not in (None, expected_mate):
                raise ValueError(
                    f"{path}: FASTQ header identifies R{observed_mate}, "
                    f"but R{expected_mate} was expected at record {record_number}"
                )

            yield FastqRecord(
                read_id=normalize_read_id(header),
                header=header,
                sequence=sequence,
                quality=quality,
            )
