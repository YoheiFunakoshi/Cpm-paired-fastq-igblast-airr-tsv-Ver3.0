from __future__ import annotations

import os
import inspect
from pathlib import Path
import shutil
import unittest
import uuid

from airr_igblast_paired.fastq import FastqRecord
from airr_igblast_paired.prepare import (
    ReadTransform,
    ensure_distinct_paths,
    prepare_paired_fastq_to_fasta,
    reverse_complement,
    transform_sequence,
)
from airr_igblast_paired.umi import CPM_R2_ANCHOR, extract_cpm_r2_umi


def record(read_id: str, sequence: str, quality: str | None = None) -> FastqRecord:
    if quality is None:
        quality = "I" * len(sequence)
    return FastqRecord(read_id=read_id, header=f"@{read_id}", sequence=sequence, quality=quality)


def cpm_r2_sequence(umi: str, insert: str = "GATTACA") -> str:
    return CPM_R2_ANCHOR + umi[:4] + "T" + umi[4:8] + "T" + umi[8:] + "TCTT" + insert


class PrepareTests(unittest.TestCase):
    def test_prepare_rejects_query_template_incompatible_with_pair_summary(self) -> None:
        root = Path(f"test_prepare_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            output = root / "queries.fasta"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "pipe-delimited R1/R2"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}_{read}",
                )
            with self.assertRaisesRegex(ValueError, r"include the \{read_id\} placeholder"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="constant|{read}",
                )
            with self.assertRaisesRegex(ValueError, r"include the \{umi\} placeholder"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}|{read}",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, r"emit UMI=\{umi\}"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}|{read}|UMI=ACGTACGTACGT|{umi}",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, r"UMI=\{umi\} after the \{read_id\}"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="UMI={umi}|{read_id}|{read}",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, "exactly one pipe-delimited UMI component"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}|{read}|UMI={umi}|UMI=FIXED",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, r"after \{read_id\}"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read}|{read_id}|UMI={umi}",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, "literal pipe-delimited R1/R2"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}|{read}|UMI={umi}|R1",
                    umi_mode="cpm-r2",
                )
            with self.assertRaisesRegex(ValueError, "format specifications"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="PAIR_ID|{read_id:.0}|{read}|UMI={umi}",
                    umi_mode="cpm-r2",
                )

            with self.assertRaisesRegex(ValueError, r"include the \{read_id\} placeholder"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="PAIR_ID|{read}",
                )
            with self.assertRaisesRegex(ValueError, r"include the \{umi\} placeholder"):
                prepare_paired_fastq_to_fasta(
                    r1,
                    r2,
                    output,
                    query_name_template="{read_id}|{read}|UMI=ACGTACGTACGT",
                    umi_mode="cpm-r2",
                )

            self.assertFalse(output.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_prepare_rejects_reserved_pipe_in_fastq_read_id(self) -> None:
        root = Path(f"test_prepare_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            output = root / "queries.fasta"
            r1.write_text("@foo|R2/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@foo|R2/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "reserved '\\|' delimiter"):
                prepare_paired_fastq_to_fasta(r1, r2, output)

            self.assertFalse(output.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("ACGTNry"), "ryNACGT")

    def test_transform_trims_before_reverse_complement(self) -> None:
        transformed = transform_sequence(
            record("read1", "AACCGGTT"),
            ReadTransform("reverse-complement", trim_left=2, trim_right=2),
        )

        self.assertEqual(transformed, "CCGG")

    def test_prepare_writes_r1_and_r2_as_separate_queries(self) -> None:
        r1_path = Path("test_prepare_R1.fastq")
        r2_path = Path("test_prepare_R2.fastq")
        out_path = Path("test_prepare.fasta")
        try:
            r1_path.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_path.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)
            output = out_path.read_text(encoding="utf-8")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.total_pairs, 1)
        self.assertEqual(stats.records_written, 2)
        self.assertIn(">read1|R1\nAACCGG\n", output)
        self.assertIn(">read1|R2\nAACCTT\n", output)

    def test_ensure_distinct_paths_rejects_relative_and_absolute_aliases(self) -> None:
        path = Path("test_prepare_path_alias.fastq")

        with self.assertRaisesRegex(ValueError, "relative_path and absolute_path"):
            ensure_distinct_paths(relative_path=path, absolute_path=path.resolve())

    @unittest.skipUnless(os.name == "nt", "Windows path comparison")
    def test_ensure_distinct_paths_is_case_insensitive_on_windows(self) -> None:
        path = Path("test_prepare_CaseAlias.fastq")
        try:
            path.write_text("test\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "original and case_alias"):
                ensure_distinct_paths(original=path, case_alias=path.with_name(path.name.lower()))
        finally:
            path.unlink(missing_ok=True)

    def test_prepare_rejects_output_equal_to_input_without_modifying_fastq(self) -> None:
        r1_path = Path("test_collision_R1.fastq")
        r2_path = Path("test_collision_R2.fastq")
        r1_text = "@read1/1\nAAAA\n+\nIIII\n"
        r2_text = "@read1/2\nCCCC\n+\nIIII\n"
        try:
            for output_path, expected_names in (
                (r1_path, "r1_path and fasta_path"),
                (r2_path, "r2_path and fasta_path"),
            ):
                with self.subTest(output_path=output_path):
                    r1_path.write_text(r1_text, encoding="utf-8")
                    r2_path.write_text(r2_text, encoding="utf-8")

                    with self.assertRaisesRegex(ValueError, expected_names):
                        prepare_paired_fastq_to_fasta(r1_path, r2_path, output_path)

                    self.assertEqual(r1_path.read_text(encoding="utf-8"), r1_text)
                    self.assertEqual(r2_path.read_text(encoding="utf-8"), r2_text)
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)

    def test_prepare_rejects_same_r1_r2_path_before_writing_output(self) -> None:
        reads_path = Path("test_same_input.fastq")
        out_path = Path("test_same_input.fasta")
        reads_text = "@read1\nAAAA\n+\nIIII\n"
        try:
            reads_path.write_text(reads_text, encoding="utf-8")
            out_path.write_text("existing output\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "r1_path and r2_path"):
                prepare_paired_fastq_to_fasta(reads_path, reads_path, out_path)

            self.assertEqual(reads_path.read_text(encoding="utf-8"), reads_text)
            self.assertEqual(out_path.read_text(encoding="utf-8"), "existing output\n")
        finally:
            reads_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

    def test_prepare_preserves_existing_output_after_midstream_failure(self) -> None:
        r1_path = Path("test_atomic_failure_R1.fastq")
        r2_path = Path("test_atomic_failure_R2.fastq")
        out_path = Path("test_atomic_failure.fasta")
        try:
            r1_path.write_text(
                "@read1/1\nAAAA\n+\nIIII\n"
                "@read2/1\nCCCC\n+\nIIII\n",
                encoding="utf-8",
            )
            r2_path.write_text(
                "@read1/2\nTTTT\n+\nIIII\n"
                "@different/2\nGGGG\n+\nIIII\n",
                encoding="utf-8",
            )
            out_path.write_text("existing output\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "read ID mismatch"):
                prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)

            self.assertEqual(out_path.read_text(encoding="utf-8"), "existing output\n")
            self.assertEqual(list(out_path.parent.glob(f".{out_path.name}.*.tmp")), [])
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

    def test_prepare_atomically_replaces_existing_output_after_success(self) -> None:
        r1_path = Path("test_atomic_success_R1.fastq")
        r2_path = Path("test_atomic_success_R2.fastq")
        out_path = Path("test_atomic_success.fasta")
        try:
            r1_path.write_text("@read1/1\nAAAA\n+\nIIII\n", encoding="utf-8")
            r2_path.write_text("@read1/2\nCCCC\n+\nIIII\n", encoding="utf-8")
            out_path.write_text("existing output\n", encoding="utf-8")

            prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)

            self.assertEqual(out_path.read_text(encoding="utf-8"), ">read1|R1\nAAAA\n>read1|R2\nGGGG\n")
            self.assertEqual(list(out_path.parent.glob(f".{out_path.name}.*.tmp")), [])
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

    def test_prepare_rejects_swapped_mate_headers(self) -> None:
        r1_path = Path("test_swapped_R1.fastq")
        r2_path = Path("test_swapped_R2.fastq")
        out_path = Path("test_swapped.fasta")
        cases = (
            ("@read1/2", "@read1/1"),
            ("@instrument:run 2:N:0:ACGT", "@instrument:run 1:N:0:ACGT"),
        )
        try:
            for r1_header, r2_header in cases:
                with self.subTest(r1_header=r1_header):
                    r1_path.write_text(f"{r1_header}\nAAAA\n+\nIIII\n", encoding="utf-8")
                    r2_path.write_text(f"{r2_header}\nCCCC\n+\nIIII\n", encoding="utf-8")
                    out_path.write_text("existing output\n", encoding="utf-8")

                    with self.assertRaisesRegex(ValueError, r"identifies R2, but R1 was expected"):
                        prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)

                    self.assertEqual(out_path.read_text(encoding="utf-8"), "existing output\n")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

    def test_prepare_allows_headers_without_mate_metadata(self) -> None:
        r1_path = Path("test_no_mate_R1.fastq")
        r2_path = Path("test_no_mate_R2.fastq")
        out_path = Path("test_no_mate.fasta")
        try:
            r1_path.write_text("@read1\nAAAA\n+\nIIII\n", encoding="utf-8")
            r2_path.write_text("@read1\nCCCC\n+\nIIII\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.records_written, 2)

    def test_prepare_can_filter_by_read_length(self) -> None:
        r1_path = Path("test_filter_R1.fastq")
        r2_path = Path("test_filter_R2.fastq")
        out_path = Path("test_filter.fasta")
        try:
            r1_path.write_text("@read1/1\nAAAA\n+\nIIII\n", encoding="utf-8")
            r2_path.write_text("@read1/2\nCCCC\n+\nIIII\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path, min_length=5)
            output = out_path.read_text(encoding="utf-8")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.records_written, 0)
        self.assertEqual(stats.skipped_too_short, 2)
        self.assertEqual(output, "")

    def test_extract_cpm_r2_umi(self) -> None:
        sequence = CPM_R2_ANCHOR + "ACGTTTGCATGATTTCTTAAA"

        self.assertEqual(extract_cpm_r2_umi(sequence), "ACGTTGCAGATT")

    def test_prepare_can_add_cpm_umi_to_query_name(self) -> None:
        r1_path = Path("test_umi_R1.fastq")
        r2_path = Path("test_umi_R2.fastq")
        out_path = Path("test_umi.fasta")
        r2_sequence = CPM_R2_ANCHOR + "ACGTTTGCATGATTTCTTAAA"
        try:
            r1_path.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_path.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(
                r1_path,
                r2_path,
                out_path,
                read_selection="r2",
                r2_transform=ReadTransform("forward"),
                query_name_template="{read_id}|{read}|UMI={umi}",
                umi_mode="cpm-r2",
            )
            output = out_path.read_text(encoding="utf-8")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.umi_extracted_pairs, 1)
        self.assertEqual(stats.umi_missing_pairs, 0)
        self.assertIn(">read1|R2|UMI=ACGTTGCAGATT\n", output)

    def test_prepare_api_has_no_pre_igblast_umi_collapse_options(self) -> None:
        parameters = inspect.signature(prepare_paired_fastq_to_fasta).parameters

        self.assertNotIn("umi_collapse", parameters)
        self.assertNotIn("umi_collapse_mismatches", parameters)
        self.assertNotIn("umi_sequence_distance", parameters)
        self.assertNotIn("umi_collapse_strategy", parameters)

    def test_prepare_retains_every_pair_and_preserves_exact_raw_umi(self) -> None:
        root = Path(f"test_prepare_all_pairs_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1_path = root / "sample_R1.fastq"
            r2_path = root / "sample_R2.fastq"
            out_path = root / "queries.fasta"
            r1_sequence = "AACCGGTTAACC"
            umis = ("AAAAAAAAAAAA", "AAAAAAAAAAAA", "AAAAAAAAAAAT")
            r1_path.write_text(
                "".join(
                    f"@read{index}/1\n{r1_sequence}\n+\n{'I' * len(r1_sequence)}\n"
                    for index in range(1, 4)
                ),
                encoding="utf-8",
            )
            r2_path.write_text(
                "".join(
                    f"@read{index}/2\n{sequence}\n+\n{'I' * len(sequence)}\n"
                    for index, sequence in enumerate(
                        (cpm_r2_sequence(umi) for umi in umis), start=1
                    )
                ),
                encoding="utf-8",
            )

            stats = prepare_paired_fastq_to_fasta(
                r1_path,
                r2_path,
                out_path,
                query_name_template="{read_id}|{read}|UMI={umi}",
                r2_transform=ReadTransform("forward"),
                umi_mode="cpm-r2",
            )
            output = out_path.read_text(encoding="utf-8")
        finally:
            shutil.rmtree(root, ignore_errors=True)

        self.assertEqual(stats.total_pairs, 3)
        self.assertEqual(stats.umi_extracted_pairs, 3)
        self.assertEqual(stats.records_written, 6)
        self.assertNotIn("umi_collapsed_pairs", vars(stats))
        self.assertNotIn("umi_duplicate_pairs_skipped", vars(stats))
        self.assertIn(">read1|R1|UMI=AAAAAAAAAAAA\n", output)
        self.assertIn(">read2|R1|UMI=AAAAAAAAAAAA\n", output)
        self.assertIn(">read3|R1|UMI=AAAAAAAAAAAT\n", output)

    def test_ambiguous_raw_umi_is_retained_but_counted_as_missing(self) -> None:
        anchor = "TATCAACGCAGAGTGGCCAT"
        r1_path = Path("ambiguous_umi_R1.fastq")
        r2_path = Path("ambiguous_umi_R2.fastq")
        out_path = Path("ambiguous_umi.fasta")
        self.addCleanup(lambda: r1_path.unlink(missing_ok=True))
        self.addCleanup(lambda: r2_path.unlink(missing_ok=True))
        self.addCleanup(lambda: out_path.unlink(missing_ok=True))
        r1_sequence = "ACGTACGT"
        r2_sequence = anchor + "AAAA" + "T" + "NAAA" + "T" + "AAAA"
        r1_path.write_text(
            f"@pair1 1:N:0:1\n{r1_sequence}\n+\n{'I' * len(r1_sequence)}\n",
            encoding="utf-8",
        )
        r2_path.write_text(
            f"@pair1 2:N:0:1\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n",
            encoding="utf-8",
        )

        stats = prepare_paired_fastq_to_fasta(
            r1_path,
            r2_path,
            out_path,
            umi_mode="cpm-r2",
            query_name_template="{read_id}|{read}|UMI={umi}",
        )

        self.assertEqual(stats.umi_extracted_pairs, 0)
        self.assertEqual(stats.umi_missing_pairs, 1)
        self.assertIn("UMI=AAAANAAAAAAA", out_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
