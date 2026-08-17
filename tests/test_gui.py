from __future__ import annotations

from pathlib import Path
import shutil
import tkinter as tk
import unittest
import uuid
from unittest.mock import Mock, patch

from airr_igblast_paired.gui import APP_TITLE, App, _normalize_run_values, _preflight_run_errors, main


def form_values() -> dict[str, object]:
    return {
        "r1": "sample_R1.fastq",
        "r2": "sample_R2.fastq",
        "out": "sample.airr.tsv",
        "query_fasta": "",
        "igblastn": "igblastn",
        "germline_db_v": "v",
        "germline_db_d": "d",
        "germline_db_j": "j",
        "auxiliary_data": "human_gl.aux",
        "organism": "human",
        "domain_system": "imgt",
        "ig_seqtype": "Ig",
        "num_threads": 4,
        "igblast_batch_size": 10000,
        "read_selection": "both",
        "r1_orientation": "forward",
        "r2_orientation": "reverse-complement",
        "trim_left_r1": 0,
        "trim_right_r1": 0,
        "trim_left_r2": 0,
        "trim_right_r2": 0,
        "min_length": 0,
        "max_n_rate": 1.0,
        "query_name_template": "{read_id}|{read}|UMI={umi}",
        "umi_anchor_max_mismatches": 2,
        "strict_ids": True,
    }


class GuiTests(unittest.TestCase):
    def test_normalize_run_values_reports_invalid_numeric_field(self) -> None:
        values = form_values()
        values["num_threads"] = "not-an-integer"

        with self.assertRaisesRegex(ValueError, "Threads には整数"):
            _normalize_run_values(values)

    def test_normalize_run_values_keeps_ver3_exact_umi_inputs_without_legacy_controls(self) -> None:
        values = _normalize_run_values(form_values())

        self.assertEqual(values["umi_anchor_max_mismatches"], 2)
        self.assertNotIn("umi_output_mode", values)
        self.assertNotIn("umi_collapse_mismatches", values)
        self.assertNotIn("umi_sequence_distance", values)

    def test_normalize_run_values_reports_out_of_range_values(self) -> None:
        values = form_values()
        values["max_n_rate"] = 1.1
        with self.assertRaisesRegex(ValueError, "0以上1以下"):
            _normalize_run_values(values)

        values = form_values()
        values["num_threads"] = 0
        with self.assertRaisesRegex(ValueError, "1以上"):
            _normalize_run_values(values)

    def test_preflight_accepts_files_db_components_and_path_igblastn(self) -> None:
        root = Path(f"test_gui_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            values = form_values()
            for key, filename in (
                ("r1", "sample_R1.fastq"),
                ("r2", "sample_R2.fastq"),
                ("igblastn", "igblastn.exe"),
                ("auxiliary_data", "human_gl.aux"),
            ):
                path = root / filename
                path.write_text("test\n", encoding="utf-8")
                values[key] = str(path)

            for key, filename in (
                ("germline_db_v", "v"),
                ("germline_db_d", "d"),
                ("germline_db_j", "j"),
            ):
                prefix = root / filename
                Path(str(prefix) + ".00.nsq").write_text("db\n", encoding="utf-8")
                values[key] = str(prefix)

            planned = (root / "result.airr.tsv",)
            with patch("airr_igblast_paired.gui.shutil.which", return_value=None):
                self.assertEqual(_preflight_run_errors(values, planned), [])

            values["igblastn"] = "igblastn"
            with patch(
                "airr_igblast_paired.gui.shutil.which",
                return_value=str(root / "igblastn.exe"),
            ):
                self.assertEqual(_preflight_run_errors(values, planned), [])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_preflight_lists_missing_files_db_and_path_collisions_in_japanese(self) -> None:
        root = Path(f"test_gui_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "same.fastq"
            r1.write_text("test\n", encoding="utf-8")
            values = form_values()
            values.update(
                {
                    "r1": str(r1),
                    "r2": str(r1),
                    "igblastn": "missing-igblastn",
                    "germline_db_v": str(root / "missing_v"),
                    "germline_db_d": str(root / "missing_d"),
                    "germline_db_j": str(root / "missing_j"),
                    "auxiliary_data": str(root / "missing.aux"),
                }
            )

            with patch("airr_igblast_paired.gui.shutil.which", return_value=None):
                errors = _preflight_run_errors(values, (r1,))

            message = "\n".join(errors)
            self.assertIn("パスが重複", message)
            self.assertIn("igblastn 実行ファイルが見つかりません", message)
            self.assertIn("V DB prefix のBLAST DB構成ファイル", message)
            self.assertIn("D DB prefix のBLAST DB構成ファイル", message)
            self.assertIn("J DB prefix のBLAST DB構成ファイル", message)
            self.assertIn("Aux file がファイルとして見つかりません", message)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_start_run_handles_tk_variable_read_error_without_starting(self) -> None:
        app = App.__new__(App)
        variable = Mock()
        variable.get.side_effect = tk.TclError("invalid integer")
        app.vars = {"num_threads": variable}
        app.running = False

        with (
            patch("airr_igblast_paired.gui.messagebox.showerror") as showerror,
            patch("airr_igblast_paired.gui.threading.Thread") as thread_class,
        ):
            app._start_run()

        showerror.assert_called_once()
        thread_class.assert_not_called()
        self.assertFalse(app.running)

    def test_start_run_handles_numeric_conversion_error_without_starting(self) -> None:
        app = App.__new__(App)
        values = form_values()
        values["num_threads"] = "invalid"
        app.vars = {key: Mock(get=Mock(return_value=value)) for key, value in values.items()}
        app.running = False

        with (
            patch("airr_igblast_paired.gui.messagebox.showerror") as showerror,
            patch("airr_igblast_paired.gui.threading.Thread") as thread_class,
        ):
            app._start_run()

        self.assertIn("Threads には整数", showerror.call_args.args[1])
        thread_class.assert_not_called()
        self.assertFalse(app.running)

    def test_start_run_restores_state_when_thread_start_fails(self) -> None:
        app = App.__new__(App)
        values = form_values()
        app.vars = {key: Mock(get=Mock(return_value=value)) for key, value in values.items()}
        app.running = False
        app.run_button = Mock()
        app._log = Mock()
        worker = Mock()
        worker.start.side_effect = RuntimeError("cannot start thread")
        planned = (Path(f"test_gui_output_{uuid.uuid4().hex}.tsv"),)

        with (
            patch(
                "airr_igblast_paired.gui.planned_cpm_output_paths",
                return_value=planned,
            ) as planned_paths,
            patch("airr_igblast_paired.gui._preflight_run_errors", return_value=[]),
            patch("airr_igblast_paired.gui.threading.Thread", return_value=worker),
            patch("airr_igblast_paired.gui.messagebox.showerror") as showerror,
        ):
            app._start_run()

        self.assertFalse(app.running)
        app.run_button.configure.assert_any_call(state="disabled")
        app.run_button.configure.assert_any_call(state="normal")
        showerror.assert_called_once()
        planned_paths.assert_called_once_with(values["out"], None)

    def test_run_pipeline_uses_single_ver3_api_without_legacy_mode_arguments(self) -> None:
        app = App.__new__(App)
        app.messages = Mock()
        multi_result = Mock(runs=(), manifest_path=Path("sample.run.json"))

        with patch(
            "airr_igblast_paired.gui.run_cpm_umi_igblast_outputs",
            return_value=multi_result,
        ) as run_pipeline:
            app._run_pipeline(form_values(), False)

        kwargs = run_pipeline.call_args.kwargs
        self.assertEqual(kwargs["umi_anchor_max_mismatches"], 2)
        self.assertNotIn("umi_output_mode", kwargs)
        self.assertNotIn("umi_collapse_mismatches", kwargs)
        self.assertNotIn("umi_sequence_distance", kwargs)
        app.messages.put.assert_called_once()

    def test_main_uses_ver3_window_title(self) -> None:
        root = Mock()
        with (
            patch("airr_igblast_paired.gui.tk.Tk", return_value=root),
            patch("airr_igblast_paired.gui.App") as app_class,
        ):
            main()

        self.assertEqual(APP_TITLE, "CPM Paired Fastq IgBLAST AIRR tsv Ver3.0")
        root.title.assert_called_once_with(APP_TITLE)
        app_class.assert_called_once_with(root)
        root.mainloop.assert_called_once_with()

    def test_close_is_refused_while_analysis_is_running(self) -> None:
        app = App.__new__(App)
        app.running = True
        app.master = Mock()

        with patch("airr_igblast_paired.gui.messagebox.showwarning") as warning:
            app._on_close()

        warning.assert_called_once()
        app.master.destroy.assert_not_called()

    def test_close_is_allowed_after_analysis_finishes(self) -> None:
        app = App.__new__(App)
        app.running = False
        app.master = Mock()

        app._on_close()

        app.master.destroy.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
