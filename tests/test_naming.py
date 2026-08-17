from __future__ import annotations

from pathlib import Path
import unittest

from airr_igblast_paired.naming import (
    PROJECT_FOLDER_NAME,
    RESULTS_FOLDER_NAME,
    default_output_tsv_path,
    fastq_stem,
    sample_name_from_fastqs,
)


class NamingTests(unittest.TestCase):
    def test_ver3_project_and_results_folder_names(self) -> None:
        self.assertEqual(PROJECT_FOLDER_NAME, "CPM Paired Fastq IgBLAST AIRR tsv Ver3.0")
        self.assertEqual(
            RESULTS_FOLDER_NAME,
            "Results of CPM Paired Fastq IgBLAST AIRR tsv Ver3.0",
        )

    def test_fastq_stem_removes_gzip_fastq_suffix(self) -> None:
        self.assertEqual(fastq_stem("sample_R1.fastq.gz"), "sample_R1")

    def test_sample_name_from_illumina_pair(self) -> None:
        sample = sample_name_from_fastqs(
            "SYNTH001_S01_L001_R1_001.fastq",
            "SYNTH001_S01_L001_R2_001.fastq",
        )

        self.assertEqual(sample, "SYNTH001_S01_L001")

    def test_sample_name_from_simple_pair(self) -> None:
        sample = sample_name_from_fastqs("synthetic-A_R1.fq.gz", "synthetic-A_R2.fq.gz")

        self.assertEqual(sample, "synthetic-A")

    def test_default_output_uses_results_folder_and_sample_name(self) -> None:
        output = default_output_tsv_path(
            "SYNTH001_S01_L001_R1_001.fastq",
            "SYNTH001_S01_L001_R2_001.fastq",
            Path("work"),
        )

        self.assertEqual(
            output,
            Path("work") / RESULTS_FOLDER_NAME / "SYNTH001_S01_L001.airr.tsv",
        )

    def test_default_output_can_add_nonlegacy_analysis_suffix(self) -> None:
        output = default_output_tsv_path(
            "SYNTH001_S01_L001_R1_001.fastq",
            "SYNTH001_S01_L001_R2_001.fastq",
            Path("work"),
            "review",
        )

        self.assertEqual(
            output,
            Path("work") / RESULTS_FOLDER_NAME / "SYNTH001_S01_L001_review.airr.tsv",
        )


if __name__ == "__main__":
    unittest.main()
