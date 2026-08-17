from __future__ import annotations

from pathlib import Path
import gzip
import unittest

from airr_igblast_paired.fastq import mate_number_from_header, normalize_read_id, read_fastq


class FastqTests(unittest.TestCase):
    def test_read_fastq_accepts_uppercase_gzip_suffix(self) -> None:
        path = Path("test_fastq_uppercase.FASTQ.GZ")
        try:
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("@read1/1\nACGT\n+\nIIII\n")
            records = list(read_fastq(path, expected_mate=1))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].sequence, "ACGT")
        finally:
            path.unlink(missing_ok=True)

    def test_normalize_read_id(self) -> None:
        self.assertEqual(normalize_read_id("@sample/1"), "sample")
        self.assertEqual(normalize_read_id("@instrument:run 2:N:0:ACGT"), "instrument:run")

    def test_mate_number_from_header_supports_common_formats(self) -> None:
        self.assertEqual(mate_number_from_header("@sample/1"), 1)
        self.assertEqual(mate_number_from_header("@sample/2"), 2)
        self.assertEqual(mate_number_from_header("@instrument:run 1:N:0:ACGT"), 1)
        self.assertEqual(mate_number_from_header("@instrument:run 2:N:0:ACGT"), 2)
        self.assertIsNone(mate_number_from_header("@read-without-mate"))

    def test_mate_number_from_header_rejects_conflicting_markers(self) -> None:
        with self.assertRaisesRegex(ValueError, "conflicting mate markers"):
            mate_number_from_header("@sample/1 2:N:0:ACGT")

    def test_read_fastq(self) -> None:
        path = Path("test_read_fastq.fastq")
        try:
            path.write_text("@read1/1\nACGT\n+\nIIII\n", encoding="utf-8")
            records = list(read_fastq(path))
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].read_id, "read1")
        self.assertEqual(records[0].sequence, "ACGT")
        self.assertEqual(records[0].quality, "IIII")

    def test_read_fastq_validates_expected_slash_mate(self) -> None:
        path = Path("test_read_fastq_wrong_slash_mate.fastq")
        try:
            path.write_text("@read1/2\nACGT\n+\nIIII\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, r"identifies R2, but R1 was expected at record 1"):
                list(read_fastq(path, expected_mate=1))
        finally:
            path.unlink(missing_ok=True)

    def test_read_fastq_validates_expected_illumina_mate(self) -> None:
        path = Path("test_read_fastq_wrong_illumina_mate.fastq")
        try:
            path.write_text("@instrument:run 2:N:0:ACGT\nACGT\n+\nIIII\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, r"identifies R2, but R1 was expected at record 1"):
                list(read_fastq(path, expected_mate=1))
        finally:
            path.unlink(missing_ok=True)

    def test_read_fastq_allows_header_without_mate_metadata(self) -> None:
        path = Path("test_read_fastq_without_mate.fastq")
        try:
            path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")

            records = list(read_fastq(path, expected_mate=2))
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual([record.read_id for record in records], ["read1"])


if __name__ == "__main__":
    unittest.main()
