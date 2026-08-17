from __future__ import annotations

import json
import inspect
import os
from pathlib import Path
import shutil
import unittest
import uuid
from unittest.mock import patch

from airr_igblast_paired.igblast import IgBlastConfig
from airr_igblast_paired.pipeline import (
    _capture_input_snapshot,
    _output_locks,
    _publish_files,
    planned_cpm_output_paths,
    run_cpm_umi_igblast_outputs,
    run_paired_igblast,
)
from airr_igblast_paired.prepare import prepare_paired_fastq_to_fasta
from airr_igblast_paired.umi import CPM_R2_ANCHOR


def cpm_r2_sequence(umi: str, insert: str = "GATTACA") -> str:
    return CPM_R2_ANCHOR + umi[:4] + "T" + umi[4:8] + "T" + umi[8:] + "TCTT" + insert


def cpm_igblast_config() -> IgBlastConfig:
    return IgBlastConfig(
        germline_db_v="v",
        germline_db_d="d",
        germline_db_j="j",
        auxiliary_data="human_gl.aux",
    )


class PipelineTests(unittest.TestCase):
    def test_run_paired_api_has_no_pre_igblast_umi_collapse_options(self) -> None:
        parameters = inspect.signature(run_paired_igblast).parameters

        self.assertNotIn("umi_collapse", parameters)
        self.assertNotIn("umi_collapse_mismatches", parameters)
        self.assertNotIn("umi_sequence_distance", parameters)
        self.assertNotIn("umi_collapse_strategy", parameters)

    def test_cpm_requires_d_database_and_auxiliary_data(self) -> None:
        base_arguments = {
            "r1_path": "sample_R1.fastq",
            "r2_path": "sample_R2.fastq",
            "output_tsv": "sample.airr.tsv",
        }
        with self.assertRaisesRegex(ValueError, "D DB prefix is required"):
            run_cpm_umi_igblast_outputs(
                **base_arguments,
                igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
            )
        with self.assertRaisesRegex(ValueError, "auxiliary data is required"):
            run_cpm_umi_igblast_outputs(
                **base_arguments,
                igblast_config=IgBlastConfig(
                    germline_db_v="v",
                    germline_db_d="d",
                    germline_db_j="j",
                ),
            )

    def test_mid_publish_failure_rolls_back_every_existing_output(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = tuple(staged / f"output{index}.txt" for index in range(3))
            final_paths = tuple(final / f"output{index}.txt" for index in range(3))
            for index, path in enumerate(staged_paths):
                path.write_text(f"new {index}\n", encoding="utf-8")
            for index, path in enumerate(final_paths):
                path.write_text(f"old {index}\n", encoding="utf-8")

            real_replace = __import__("os").replace
            published_replaces = 0

            def fail_second_publish(source: str | Path, destination: str | Path) -> None:
                nonlocal published_replaces
                source_path = Path(source)
                if source_path.name.endswith(".partial"):
                    published_replaces += 1
                    if published_replaces == 2:
                        raise OSError("injected mid-publish failure")
                real_replace(source, destination)

            with patch("airr_igblast_paired.pipeline.os.replace", side_effect=fail_second_publish):
                with self.assertRaisesRegex(OSError, "injected mid-publish failure"):
                    _publish_files(staged_paths, final_paths, overwrite=True)

            for index, path in enumerate(final_paths):
                self.assertEqual(path.read_text(encoding="utf-8"), f"old {index}\n")
            self.assertEqual(list(final.glob(".*.partial")), [])
            self.assertEqual(list(final.glob(".*.backup")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_keyboard_interrupt_during_publish_rolls_back_every_existing_output(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = tuple(staged / f"output{index}.txt" for index in range(2))
            final_paths = tuple(final / f"output{index}.txt" for index in range(2))
            for index, path in enumerate(staged_paths):
                path.write_text(f"new {index}\n", encoding="utf-8")
            for index, path in enumerate(final_paths):
                path.write_text(f"old {index}\n", encoding="utf-8")

            real_replace = __import__("os").replace
            published_replaces = 0

            def interrupt_second_publish(source: str | Path, destination: str | Path) -> None:
                nonlocal published_replaces
                source_path = Path(source)
                if source_path.name.endswith(".partial"):
                    published_replaces += 1
                    if published_replaces == 2:
                        raise KeyboardInterrupt("injected Ctrl+C")
                real_replace(source, destination)

            with patch("airr_igblast_paired.pipeline.os.replace", side_effect=interrupt_second_publish):
                with self.assertRaisesRegex(KeyboardInterrupt, r"injected Ctrl\+C"):
                    _publish_files(staged_paths, final_paths, overwrite=True)

            for index, path in enumerate(final_paths):
                self.assertEqual(path.read_text(encoding="utf-8"), f"old {index}\n")
            self.assertEqual(list(final.glob(".*.partial")), [])
            self.assertEqual(list(final.glob(".*.backup")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_completion_marker_is_backed_up_first_and_published_last(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace
            events: list[tuple[str, str]] = []

            def record_replace(source: str | Path, destination: str | Path) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if destination_path.name.endswith(".backup"):
                    events.append(("backup", source_path.name))
                elif source_path.name.endswith(".partial"):
                    events.append(("publish", destination_path.name))
                real_replace(source, destination)

            with patch("airr_igblast_paired.pipeline.os.replace", side_effect=record_replace):
                _publish_files(
                    staged_paths,
                    final_paths,
                    overwrite=True,
                    completion_marker=final_paths[-1],
                )

            backup_events = [name for kind, name in events if kind == "backup"]
            publish_events = [name for kind, name in events if kind == "publish"]
            self.assertEqual(backup_events[0], "sample.run.json")
            self.assertEqual(publish_events[-1], "sample.run.json")
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_incomplete_rollback_never_restores_completion_marker(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace
            published_replaces = 0

            def fail_publish_and_one_restore(source: str | Path, destination: str | Path) -> None:
                nonlocal published_replaces
                source_path = Path(source)
                destination_path = Path(destination)
                if source_path.name.endswith(".partial"):
                    published_replaces += 1
                    if published_replaces == 2:
                        raise OSError("injected publish failure")
                if source_path.name.endswith(".backup") and destination_path == final_paths[0]:
                    raise OSError("injected rollback failure")
                real_replace(source, destination)

            with patch(
                "airr_igblast_paired.pipeline.os.replace",
                side_effect=fail_publish_and_one_restore,
            ):
                with self.assertRaisesRegex(RuntimeError, "rollback was incomplete"):
                    _publish_files(
                        staged_paths,
                        final_paths,
                        overwrite=True,
                        completion_marker=final_paths[-1],
                    )

            self.assertFalse(final_paths[-1].exists())
            self.assertNotEqual(final_paths[0].read_text(encoding="utf-8"), "old sample.airr.tsv\n")
            self.assertTrue(list(final.glob(".*.backup")))
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_interrupt_after_successful_backup_rename_restores_old_complete_set(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace

            def interrupt_after_backup_rename(
                source: str | Path,
                destination: str | Path,
            ) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if source_path == final_paths[0] and destination_path.name.endswith(".backup"):
                    real_replace(source, destination)
                    raise KeyboardInterrupt("injected after successful backup rename")
                real_replace(source, destination)

            with patch(
                "airr_igblast_paired.pipeline.os.replace",
                side_effect=interrupt_after_backup_rename,
            ):
                with self.assertRaisesRegex(
                    KeyboardInterrupt,
                    "injected after successful backup rename",
                ):
                    _publish_files(
                        staged_paths,
                        final_paths,
                        overwrite=True,
                        completion_marker=final_paths[-1],
                    )

            for path in final_paths:
                self.assertEqual(path.read_text(encoding="utf-8"), f"old {path.name}\n")
            self.assertEqual(list(final.glob(".*.partial")), [])
            self.assertEqual(list(final.glob(".*.backup")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_interrupt_after_successful_publish_rename_restores_old_complete_set(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace

            def interrupt_after_publish_rename(
                source: str | Path,
                destination: str | Path,
            ) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if source_path.name.endswith(".partial") and destination_path == final_paths[0]:
                    real_replace(source, destination)
                    raise KeyboardInterrupt("injected after successful publish rename")
                real_replace(source, destination)

            with patch(
                "airr_igblast_paired.pipeline.os.replace",
                side_effect=interrupt_after_publish_rename,
            ):
                with self.assertRaisesRegex(
                    KeyboardInterrupt,
                    "injected after successful publish rename",
                ):
                    _publish_files(
                        staged_paths,
                        final_paths,
                        overwrite=True,
                        completion_marker=final_paths[-1],
                    )

            for path in final_paths:
                self.assertEqual(path.read_text(encoding="utf-8"), f"old {path.name}\n")
            self.assertEqual(list(final.glob(".*.partial")), [])
            self.assertEqual(list(final.glob(".*.backup")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_interrupt_after_completion_marker_commit_keeps_new_complete_set(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace
            real_unlink = Path.unlink

            def interrupt_after_marker_commit(
                source: str | Path,
                destination: str | Path,
            ) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if source_path.name.endswith(".partial") and destination_path == final_paths[-1]:
                    real_replace(source, destination)
                    raise KeyboardInterrupt("injected after completion-marker commit")
                real_replace(source, destination)

            def deny_marker_unlink(path: Path, *args, **kwargs) -> None:
                if path == final_paths[-1]:
                    raise PermissionError("injected marker unlink denial")
                real_unlink(path, *args, **kwargs)

            with patch(
                "airr_igblast_paired.pipeline.os.replace",
                side_effect=interrupt_after_marker_commit,
            ), patch(
                "airr_igblast_paired.pipeline.Path.unlink",
                autospec=True,
                side_effect=deny_marker_unlink,
            ):
                with self.assertRaisesRegex(
                    KeyboardInterrupt,
                    "injected after completion-marker commit",
                ):
                    _publish_files(
                        staged_paths,
                        final_paths,
                        overwrite=True,
                        completion_marker=final_paths[-1],
                    )

            for path in final_paths:
                self.assertEqual(path.read_text(encoding="utf-8"), f"new {path.name}\n")
            self.assertEqual(list(final.glob(".*.partial")), [])
            self.assertEqual(list(final.glob(".*.backup")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_missing_recorded_backup_does_not_restore_completion_marker(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            staged = root / "staged"
            final = root / "final"
            staged.mkdir(parents=True)
            final.mkdir()
            staged_paths = (
                staged / "sample.airr.tsv",
                staged / "sample.integrated_counts.tsv",
                staged / "sample.run.json",
            )
            final_paths = (
                final / "sample.airr.tsv",
                final / "sample.integrated_counts.tsv",
                final / "sample.run.json",
            )
            for path in staged_paths:
                path.write_text(f"new {path.name}\n", encoding="utf-8")
            for path in final_paths:
                path.write_text(f"old {path.name}\n", encoding="utf-8")

            real_replace = __import__("os").replace

            def remove_backup_then_interrupt(
                source: str | Path,
                destination: str | Path,
            ) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if source_path == final_paths[0] and destination_path.name.endswith(".backup"):
                    real_replace(source, destination)
                    destination_path.unlink()
                    raise KeyboardInterrupt("injected recorded-backup loss")
                real_replace(source, destination)

            with patch(
                "airr_igblast_paired.pipeline.os.replace",
                side_effect=remove_backup_then_interrupt,
            ):
                with self.assertRaisesRegex(RuntimeError, "rollback was incomplete"):
                    _publish_files(
                        staged_paths,
                        final_paths,
                        overwrite=True,
                        completion_marker=final_paths[-1],
                    )

            self.assertFalse(final_paths[-1].exists())
            self.assertFalse(final_paths[0].exists())
            self.assertEqual(
                final_paths[1].read_text(encoding="utf-8"),
                "old sample.integrated_counts.tsv\n",
            )
            self.assertTrue(list(final.glob(".sample.run.json.*.backup")))
            self.assertEqual(list(final.glob(".*.partial")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_aliases_for_same_cpm_outputs_share_locks(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            first = planned_cpm_output_paths(root / "sample.airr.tsv")
            alias = planned_cpm_output_paths(root / "sample_umiSeq5.airr.tsv")
            self.assertEqual(
                {path.resolve(strict=False) for path in first},
                {path.resolve(strict=False) for path in alias},
            )
            self.assertEqual(first[-1], root / "sample.run.json")
            self.assertEqual(alias[-1], root / "sample.run.json")

            with _output_locks(first):
                with self.assertRaisesRegex(RuntimeError, "overlapping output"):
                    with _output_locks(alias):
                        self.fail("overlapping output locks must not be acquired")

            self.assertEqual(list(root.glob("*.lock")), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_query_and_airr_output_collision_is_rejected_before_reading(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            output = root / "sample.airr.tsv"
            r1_bytes = b"@read1/1\nAACCGG\n+\nIIIIII\n"
            r2_bytes = b"@read1/2\nAAGGTT\n+\nIIIIII\n"
            r1.write_bytes(r1_bytes)
            r2.write_bytes(r2_bytes)

            with self.assertRaisesRegex(ValueError, "paths must be distinct"):
                run_paired_igblast(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=output,
                    query_fasta=output,
                    igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                    work_dir=root / "work",
                )

            self.assertEqual(r1.read_bytes(), r1_bytes)
            self.assertEqual(r2.read_bytes(), r2_bytes)
            self.assertFalse(output.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_cpm_outputs_cannot_replace_igblast_reference_or_executable_files(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")
            executable = root / "igblast.exe"
            v_prefix = root / "V"
            v_component = root / "V.nsq"
            auxiliary = root / "human_gl.aux"
            executable.write_bytes(b"executable")
            v_component.write_bytes(b"database")
            auxiliary.write_bytes(b"auxiliary")

            protected_cases = (
                (executable, "IgBLAST executable"),
                (v_component, "V DB component"),
                (auxiliary, "IgBLAST auxiliary data"),
            )
            for protected_path, expected_label in protected_cases:
                before = protected_path.read_bytes()
                with self.subTest(path=protected_path.name):
                    with self.assertRaisesRegex(ValueError, expected_label):
                        run_cpm_umi_igblast_outputs(
                            r1_path=r1,
                            r2_path=r2,
                            output_tsv=root / "results" / "sample.airr.tsv",
                            query_fasta=protected_path,
                            igblast_config=IgBlastConfig(
                                igblastn=str(executable),
                                germline_db_v=str(v_prefix),
                                germline_db_d=str(root / "d"),
                                germline_db_j=str(root / "j"),
                                auxiliary_data=str(auxiliary),
                            ),
                        )
                    self.assertEqual(protected_path.read_bytes(), before)
            self.assertFalse((root / "results").exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_work_dir_stages_all_outputs_before_publishing(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            def fake_run_igblast(query_fasta: Path, output_tsv: Path, _: IgBlastConfig) -> list[str]:
                self.assertTrue(str(query_fasta).startswith(str(work.resolve())))
                self.assertTrue(str(output_tsv).startswith(str(work.resolve())))
                self.assertNotEqual(output_tsv, out)
                output_tsv.write_text(
                    "sequence_id\tsequence\tv_call\tj_call\tjunction_aa\tproductive\n"
                    "read1|R1\tAACCGG\tIGHV1\tIGHJ4\tCARYW\tT\n"
                    "read1|R2\tAACCTT\tIGHV1\tIGHJ4\tCARYW\tT\n",
                    encoding="utf-8",
                )
                return ["igblastn", "-query", str(query_fasta), "-out", str(output_tsv)]

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fake_run_igblast):
                result = run_paired_igblast(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=out,
                    query_fasta=query,
                    igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                    work_dir=work,
                )

            self.assertEqual(result.output_tsv, out)
            self.assertEqual(result.query_fasta, query)
            self.assertTrue(out.exists())
            self.assertTrue(query.exists())
            self.assertTrue((root / "results" / "sample.R1.airr.tsv").exists())
            self.assertTrue((root / "results" / "sample.R2.airr.tsv").exists())
            self.assertTrue((root / "results" / "sample.integrated.tsv").exists())
            self.assertTrue((root / "results" / "sample.integrated_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.integrated_counts.xlsx").exists())
            self.assertTrue((root / "results" / "sample.final_productive_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.final_productive_counts.xlsx").exists())
            self.assertTrue((root / "results" / "sample.umi_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.umi_counts.xlsx").exists())
            self.assertTrue((root / "results" / "sample.final_productive_umi_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.final_productive_umi_counts.xlsx").exists())
            self.assertEqual(result.pair_summary_stats.total_pairs, 1)
            self.assertEqual(result.pair_summary_stats.included_in_counts, 1)
            self.assertEqual(result.pair_summary_stats.unique_final_clonotypes, 1)
            self.assertEqual(result.pair_summary_stats.final_productive_included_in_counts, 1)
            self.assertEqual(result.pair_summary_stats.unique_final_productive_clonotypes, 1)
            self.assertIn(">read1|R2\nAACCTT\n", query.read_text(encoding="utf-8"))
            self.assertEqual(list(work.iterdir()), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_work_dir_cleans_partial_outputs_after_failure(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            def fake_run_igblast(_: Path, output_tsv: Path, __: IgBlastConfig) -> list[str]:
                output_tsv.write_text("partial\n", encoding="utf-8")
                raise RuntimeError("IgBLAST failed")

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fake_run_igblast):
                with self.assertRaisesRegex(RuntimeError, "IgBLAST failed"):
                    run_paired_igblast(
                        r1_path=r1,
                        r2_path=r2,
                        output_tsv=out,
                        query_fasta=query,
                        igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                        work_dir=work,
                    )

            self.assertFalse(out.exists())
            self.assertFalse(query.exists())
            self.assertEqual(list(work.iterdir()), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_nonempty_query_with_missing_airr_data_rows_is_not_published(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            for airr_text in (
                "\n",
                "sequence_id\tv_call\tj_call\tjunction_aa\tproductive\n",
                "sequence_id\tv_call\tj_call\tjunction_aa\tproductive\n"
                "read1|R1\tIGHV1\tIGHJ4\tCARYW\tT\n",
            ):
                with self.subTest(airr_text=repr(airr_text)):
                    shutil.rmtree(root / "results", ignore_errors=True)

                    def fake_run_igblast(
                        _query_fasta: Path,
                        output_tsv: Path,
                        _config: IgBlastConfig,
                    ) -> list[str]:
                        output_tsv.write_text(airr_text, encoding="utf-8")
                        return ["igblastn"]

                    with patch(
                        "airr_igblast_paired.pipeline.run_igblast",
                        side_effect=fake_run_igblast,
                    ):
                        with self.assertRaisesRegex(RuntimeError, "row count does not match"):
                            run_paired_igblast(
                                r1_path=r1,
                                r2_path=r2,
                                output_tsv=out,
                                igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                                work_dir=root / "work",
                            )

                    self.assertFalse(out.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_planned_cpm_outputs_are_single_base_analysis_with_umi_tables(self) -> None:
        root = Path("planned")
        paths = planned_cpm_output_paths(
            root / "sample_umiSeq5.airr.tsv",
            root / "sample_umiNoCollapse.queries.fasta",
        )

        self.assertEqual(
            [path.name for path in paths],
            [
                "sample.airr.tsv",
                "sample.R1.airr.tsv",
                "sample.R2.airr.tsv",
                "sample.integrated.tsv",
                "sample.integrated_counts.tsv",
                "sample.integrated_counts.xlsx",
                "sample.final_productive_counts.tsv",
                "sample.final_productive_counts.xlsx",
                "sample.umi_counts.tsv",
                "sample.umi_counts.xlsx",
                "sample.final_productive_umi_counts.tsv",
                "sample.final_productive_umi_counts.xlsx",
                "sample.queries.fasta",
                "sample.run.json",
            ],
        )
        self.assertEqual(len(paths), len(set(paths)))

    def test_cpm_runs_prepare_and_igblast_once_without_legacy_suffixes(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            seen_outputs: list[Path] = []

            def fake_run_igblast(query_fasta: Path, output_tsv: Path, _: IgBlastConfig) -> list[str]:
                seen_outputs.append(output_tsv)
                self.assertIn("UMI=ACGTACGTACGT", query_fasta.read_text(encoding="utf-8"))
                output_tsv.write_text(
                    "sequence_id\tsequence\tv_call\tj_call\tjunction_aa\tproductive\n"
                    "read1|R1|UMI=ACGTACGTACGT\tAACCGG\tIGHV1\tIGHJ4\tCARYW\tT\n"
                    "read1|R2|UMI=ACGTACGTACGT\tAACCTT\tIGHV1\tIGHJ4\tCARYW\tT\n",
                    encoding="utf-8",
                )
                return ["igblastn", "-query", str(query_fasta), "-out", str(output_tsv)]

            with patch(
                "airr_igblast_paired.pipeline.prepare_paired_fastq_to_fasta",
                wraps=prepare_paired_fastq_to_fasta,
            ) as prepare_mock, patch(
                "airr_igblast_paired.pipeline.run_igblast",
                side_effect=fake_run_igblast,
            ):
                result = run_cpm_umi_igblast_outputs(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=out,
                    query_fasta=query,
                    igblast_config=cpm_igblast_config(),
                    work_dir=root / "work",
                )

            self.assertEqual(prepare_mock.call_count, 1)
            self.assertEqual([run.analysis_suffix for run in result.runs], [""])
            self.assertEqual(len(seen_outputs), 1)
            self.assertTrue(all(path.parent != root / "results" for path in seen_outputs))
            self.assertTrue(all(path.name == "analysis.airr.tsv" for path in seen_outputs))
            self.assertEqual(result.runs[0].result.output_tsv, out)
            self.assertTrue((root / "results" / "sample.queries.fasta").exists())
            self.assertTrue((root / "results" / "sample.integrated_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.final_productive_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.umi_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.umi_counts.xlsx").exists())
            self.assertTrue((root / "results" / "sample.final_productive_umi_counts.tsv").exists())
            self.assertTrue((root / "results" / "sample.final_productive_umi_counts.xlsx").exists())
            published_names = [path.name for path in (root / "results").iterdir()]
            self.assertFalse(any("umiSeq" in name or "umiNoCollapse" in name for name in published_names))
            self.assertEqual(result.runs[0].result.stats.umi_extracted_pairs, 1)
            self.assertNotIn("umi_duplicate_pairs_skipped", vars(result.runs[0].result.stats))

            self.assertEqual(result.manifest_path, root / "results" / "sample.run.json")
            manifest_text = result.manifest_path.read_text(encoding="utf-8")
            manifest = json.loads(manifest_text)
            self.assertEqual(manifest["manifest_schema_version"], 2)
            self.assertEqual(
                manifest["counting_semantics"],
                "cpm_v3_read_pair_and_exact_raw_umi_per_clonotype_v1",
            )
            self.assertTrue(manifest["completed_utc"].endswith("Z"))
            self.assertEqual(manifest["manifest_path"], str(result.manifest_path.resolve()))
            self.assertEqual(manifest["inputs"]["r1"]["path"], str(r1.resolve()))
            self.assertEqual(manifest["inputs"]["r1"]["size"], r1.stat().st_size)
            self.assertEqual(set(manifest["inputs"]["r1"]), {"path", "size", "mtime"})
            self.assertEqual(manifest["inputs"]["r2"]["path"], str(r2.resolve()))
            self.assertEqual(set(manifest["inputs"]["r2"]), {"path", "size", "mtime"})
            self.assertNotIn("sha256", manifest_text.lower())
            self.assertEqual(manifest["settings"]["umi"]["mode"], "cpm-r2")
            self.assertEqual(
                manifest["settings"]["umi"]["counting_unit"],
                "exact_raw_umi_per_clonotype",
            )
            self.assertNotIn("collapse", manifest["settings"]["umi"])
            self.assertNotIn("output_mode", manifest["settings"]["umi"])
            self.assertNotIn("sequence_distance", manifest["settings"]["umi"])
            self.assertNotIn("collapse_mismatches", manifest["settings"]["umi"])
            self.assertEqual(manifest["settings"]["qc"]["min_length"], 0)
            self.assertEqual(manifest["settings"]["igblast"]["germline_db_v"], "v")
            self.assertEqual(
                [mode["analysis_suffix"] for mode in manifest["modes"]],
                [""],
            )
            self.assertNotIn("umi_collapse", manifest["modes"][0])
            self.assertEqual(manifest["modes"][0]["umi_counting"], "exact_raw_umi_per_clonotype")
            self.assertEqual(
                {output["name"] for output in manifest["modes"][0]["outputs"]},
                {
                    "output_tsv",
                    "r1_tsv",
                    "r2_tsv",
                    "integrated_tsv",
                    "counts_tsv",
                    "counts_xlsx",
                    "final_productive_counts_tsv",
                    "final_productive_counts_xlsx",
                    "umi_counts_tsv",
                    "umi_counts_xlsx",
                    "final_productive_umi_counts_tsv",
                    "final_productive_umi_counts_xlsx",
                    "query_fasta",
                },
            )
            for mode in manifest["modes"]:
                self.assertTrue(mode["label"])
                self.assertIn("prepare", mode["stats"])
                self.assertIn("pair_summary", mode["stats"])
                self.assertTrue(mode["outputs"])
                for output in mode["outputs"]:
                    output_path = Path(output["path"])
                    self.assertTrue(output_path.is_file())
                    self.assertEqual(output["size"], output_path.stat().st_size)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_existing_outputs_require_explicit_overwrite(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            existing = root / "results" / "sample.airr.tsv"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            existing.parent.mkdir()
            existing.write_text("old result\n", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                run_cpm_umi_igblast_outputs(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=out,
                    igblast_config=cpm_igblast_config(),
                )

            self.assertEqual(existing.read_text(encoding="utf-8"), "old result\n")
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_existing_completion_manifest_requires_explicit_overwrite(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            manifest = root / "results" / "sample.run.json"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            manifest.parent.mkdir()
            manifest.write_text('{"old": true}\n', encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                run_cpm_umi_igblast_outputs(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=out,
                    igblast_config=cpm_igblast_config(),
                )

            self.assertEqual(manifest.read_text(encoding="utf-8"), '{"old": true}\n')
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_cpm_failure_publishes_nothing(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            def fake_run_igblast(_: Path, output_tsv: Path, __: IgBlastConfig) -> list[str]:
                output_tsv.write_text("partial\n", encoding="utf-8")
                raise RuntimeError("single analysis failed")

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fake_run_igblast):
                with self.assertRaisesRegex(RuntimeError, "single analysis failed"):
                    run_cpm_umi_igblast_outputs(
                        r1_path=r1,
                        r2_path=r2,
                        output_tsv=out,
                        query_fasta=query,
                        igblast_config=cpm_igblast_config(),
                        work_dir=work,
                    )

            self.assertEqual(list((root / "results").iterdir()), [])
            self.assertFalse((root / "results" / "sample.run.json").exists())
            self.assertEqual(list(work.iterdir()), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_input_replacement_before_publish_publishes_nothing(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            initial_size = r1.stat().st_size
            initial_mtime_ns = r1.stat().st_mtime_ns
            igblast_calls = 0

            def fake_run_igblast(
                _query_fasta: Path,
                output_tsv: Path,
                _config: IgBlastConfig,
            ) -> list[str]:
                nonlocal igblast_calls
                igblast_calls += 1
                output_tsv.write_text(
                    "sequence_id\tsequence\tv_call\tj_call\tjunction_aa\tproductive\n"
                    "read1|R1|UMI=ACGTACGTACGT\tAACCGG\tIGHV1\tIGHJ4\tCARYW\tT\n"
                    "read1|R2|UMI=ACGTACGTACGT\tAACCTT\tIGHV1\tIGHJ4\tCARYW\tT\n",
                    encoding="utf-8",
                )
                return ["igblastn"]

            def mutate_before_publish(message: str) -> None:
                if message != "Publishing completed CPM outputs...":
                    return
                replacement = root / "replacement_R1.fastq"
                replacement.write_text("@read1/1\nTTTTTT\n+\nIIIIII\n", encoding="utf-8")
                self.assertEqual(replacement.stat().st_size, initial_size)
                os.utime(
                    replacement,
                    ns=(replacement.stat().st_atime_ns, initial_mtime_ns + 1_000_000_000),
                )
                os.replace(replacement, r1)

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fake_run_igblast):
                with self.assertRaisesRegex(RuntimeError, "R1 FASTQ changed during analysis"):
                    run_cpm_umi_igblast_outputs(
                        r1_path=r1,
                        r2_path=r2,
                        output_tsv=out,
                        igblast_config=cpm_igblast_config(),
                        work_dir=work,
                        progress_callback=mutate_before_publish,
                    )

            self.assertEqual(igblast_calls, 1)
            self.assertFalse((root / "results" / "sample.run.json").exists())
            self.assertEqual(list((root / "results").iterdir()), [])
            self.assertEqual(list(work.iterdir()), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_input_retargeted_to_output_during_snapshot_is_rejected(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            colliding_output = root / "results" / "sample.airr.tsv"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            colliding_output.parent.mkdir()
            colliding_output.write_text("old protected output\n", encoding="utf-8")

            capture_calls = 0

            def retarget_r1_before_capture(path: str | Path):
                nonlocal capture_calls
                capture_calls += 1
                if capture_calls == 1:
                    r1.unlink()
                    os.link(colliding_output, r1)
                return _capture_input_snapshot(path)

            with patch(
                "airr_igblast_paired.pipeline._capture_input_snapshot",
                side_effect=retarget_r1_before_capture,
            ):
                with self.assertRaisesRegex(ValueError, "paths must be distinct"):
                    run_cpm_umi_igblast_outputs(
                        r1_path=r1,
                        r2_path=r2,
                        output_tsv=out,
                        igblast_config=cpm_igblast_config(),
                        work_dir=work,
                        overwrite=True,
                    )

            self.assertEqual(
                colliding_output.read_text(encoding="utf-8"),
                "old protected output\n",
            )
            self.assertFalse((root / "results" / "sample.run.json").exists())
            self.assertFalse(work.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_overwrite_failure_preserves_existing_complete_outputs(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            old_output = root / "results" / "sample.airr.tsv"
            old_query = root / "results" / "sample.queries.fasta"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_sequence = cpm_r2_sequence("ACGTACGTACGT")
            r2.write_text(f"@read1/2\n{r2_sequence}\n+\n{'I' * len(r2_sequence)}\n", encoding="utf-8")
            old_output.parent.mkdir()
            old_output.write_text("old AIRR\n", encoding="utf-8")
            old_query.write_text("old query\n", encoding="utf-8")

            def fail_run(_: Path, output_tsv: Path, __: IgBlastConfig) -> list[str]:
                output_tsv.write_text("new partial\n", encoding="utf-8")
                raise RuntimeError("injected failure")

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fail_run):
                with self.assertRaisesRegex(RuntimeError, "injected failure"):
                    run_cpm_umi_igblast_outputs(
                        r1_path=r1,
                        r2_path=r2,
                        output_tsv=out,
                        query_fasta=query,
                        igblast_config=cpm_igblast_config(),
                        work_dir=root / "work",
                        overwrite=True,
                    )

            self.assertEqual(old_output.read_text(encoding="utf-8"), "old AIRR\n")
            self.assertEqual(old_query.read_text(encoding="utf-8"), "old query\n")
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
