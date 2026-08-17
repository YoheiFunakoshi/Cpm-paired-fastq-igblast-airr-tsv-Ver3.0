from __future__ import annotations

from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

from airr_igblast_paired import __version__
from airr_igblast_paired.cli import build_parser
from airr_igblast_paired.pipeline import MultiPipelineResult, NamedPipelineResult, PipelineResult
from airr_igblast_paired.prepare import PrepareStats


class CliTests(unittest.TestCase):
    def test_run_reports_both_exact_umi_family_views(self) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "run",
                "--r1",
                "sample_R1.fastq",
                "--r2",
                "sample_R2.fastq",
                "--germline-db-v",
                "v",
                "--germline-db-d",
                "d",
                "--germline-db-j",
                "j",
                "--auxiliary-data",
                "aux",
                "--out",
                "sample.airr.tsv",
            ]
        )
        result = PipelineResult(
            stats=PrepareStats(),
            command=["igblastn"],
            query_fasta=None,
            output_tsv=Path("sample.airr.tsv"),
            exact_umi_family_counts_tsv=Path("sample.exact_umi_family_counts.tsv"),
            exact_umi_family_counts_xlsx=Path("sample.exact_umi_family_counts.xlsx"),
            final_productive_exact_umi_family_counts_tsv=Path(
                "sample.final_productive_exact_umi_family_counts.tsv"
            ),
            final_productive_exact_umi_family_counts_xlsx=Path(
                "sample.final_productive_exact_umi_family_counts.xlsx"
            ),
        )
        multi_result = MultiPipelineResult(
            runs=(NamedPipelineResult("test", "", result),),
            manifest_path=Path("sample.run.json"),
        )
        stderr = StringIO()

        with (
            patch(
                "airr_igblast_paired.cli.run_cpm_umi_igblast_outputs",
                return_value=multi_result,
            ),
            redirect_stderr(stderr),
        ):
            self.assertEqual(args.func(args), 0)

        output = stderr.getvalue()
        self.assertIn("sample.exact_umi_family_counts.tsv", output)
        self.assertIn("sample.exact_umi_family_counts.xlsx", output)
        self.assertIn("sample.final_productive_exact_umi_family_counts.tsv", output)
        self.assertIn("sample.final_productive_exact_umi_family_counts.xlsx", output)

    def test_ver3_cli_defaults_to_one_annotation_run_with_exact_raw_umi(self) -> None:
        parser = build_parser()

        self.assertEqual(parser.prog, "cpm-paired-fastq-igblast-airr-tsv-v3")
        self.assertEqual(__version__, "3.0.1")

        prepare_args = parser.parse_args(
            [
                "prepare",
                "--r1",
                "sample_R1.fastq",
                "--r2",
                "sample_R2.fastq",
                "--out-fasta",
                "sample.queries.fasta",
            ]
        )

        self.assertEqual(prepare_args.read_selection, "both")
        self.assertEqual(prepare_args.umi_mode, "cpm-r2")
        self.assertEqual(prepare_args.umi_anchor_max_mismatches, 2)
        self.assertFalse(hasattr(prepare_args, "umi_collapse"))
        self.assertFalse(hasattr(prepare_args, "umi_collapse_mismatches"))
        self.assertFalse(hasattr(prepare_args, "umi_sequence_distance"))
        self.assertFalse(hasattr(prepare_args, "umi_collapse_strategy"))
        self.assertEqual(prepare_args.query_name_template, "{read_id}|{read}|UMI={umi}")

        run_args = parser.parse_args(
            [
                "run",
                "--r1",
                "sample_R1.fastq",
                "--r2",
                "sample_R2.fastq",
                "--germline-db-v",
                "v",
                "--germline-db-d",
                "d",
                "--germline-db-j",
                "j",
                "--auxiliary-data",
                "aux",
                "--out",
                "sample.airr.tsv",
            ]
        )
        self.assertEqual(run_args.read_selection, "both")
        self.assertEqual(run_args.umi_anchor_max_mismatches, 2)
        self.assertFalse(hasattr(run_args, "umi_output_mode"))
        self.assertFalse(hasattr(run_args, "umi_collapse_mismatches"))
        self.assertFalse(hasattr(run_args, "umi_sequence_distance"))
        self.assertEqual(run_args.igblast_batch_size, 10000)
        self.assertFalse(run_args.overwrite)

        multi_result = Mock(runs=(), manifest_path=Path("sample.run.json"))
        with (
            patch(
                "airr_igblast_paired.cli.run_cpm_umi_igblast_outputs",
                return_value=multi_result,
            ) as run_pipeline,
            redirect_stderr(StringIO()),
        ):
            self.assertEqual(run_args.func(run_args), 0)
        run_kwargs = run_pipeline.call_args.kwargs
        self.assertEqual(run_kwargs["umi_anchor_max_mismatches"], 2)
        self.assertNotIn("umi_output_mode", run_kwargs)
        self.assertNotIn("umi_collapse_mismatches", run_kwargs)
        self.assertNotIn("umi_sequence_distance", run_kwargs)

        overwrite_args = parser.parse_args(
            [
                "run",
                "--r1",
                "sample_R1.fastq",
                "--r2",
                "sample_R2.fastq",
                "--germline-db-v",
                "v",
                "--germline-db-d",
                "d",
                "--germline-db-j",
                "j",
                "--auxiliary-data",
                "aux",
                "--out",
                "sample.airr.tsv",
                "--overwrite",
            ]
        )
        self.assertTrue(overwrite_args.overwrite)

    def test_ver3_run_rejects_legacy_collapse_options(self) -> None:
        parser = build_parser()
        base_args = [
            "run",
            "--r1",
            "sample_R1.fastq",
            "--r2",
            "sample_R2.fastq",
            "--germline-db-v",
            "v",
            "--germline-db-d",
            "d",
            "--germline-db-j",
            "j",
            "--auxiliary-data",
            "aux",
            "--out",
            "sample.airr.tsv",
        ]
        for legacy_option, value in (
            ("--umi-output-mode", "both"),
            ("--umi-collapse-mismatches", "1"),
            ("--umi-sequence-distance", "5"),
        ):
            with (
                self.subTest(option=legacy_option),
                redirect_stderr(StringIO()),
                self.assertRaises(SystemExit),
            ):
                parser.parse_args([*base_args, legacy_option, value])


if __name__ == "__main__":
    unittest.main()
