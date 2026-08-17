from __future__ import annotations

import json
from pathlib import Path
import shutil
import unittest
import uuid
from types import SimpleNamespace
from unittest.mock import patch

from airr_igblast_paired.igblast import (
    IgBlastConfig,
    _copy_file_if_needed,
    _igblast_runtime_context,
    _stage_windows_igblast_config,
    _stage_windows_igblast_resources,
    _subprocess_run_options,
    build_igblast_command,
    run_igblast,
    run_igblast_batched,
    validate_extra_igblast_args,
)


class IgBlastTests(unittest.TestCase):
    def test_extra_args_cannot_override_pipeline_managed_flags(self) -> None:
        for token in ("-out", "--query", "-outfmt=7", "-germline_db_V"):
            with self.subTest(token=token):
                with self.assertRaisesRegex(ValueError, "pipeline-managed flag"):
                    validate_extra_igblast_args([token, "unsafe-value"])

        validate_extra_igblast_args(["-evalue", "1e-5"])

    def test_copy_file_updates_same_size_different_content_atomically(self) -> None:
        root = Path(f"test_igblast_copy_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            source = root / "source.bin"
            target = root / "target.bin"
            source.write_bytes(b"new!")
            target.write_bytes(b"old!")
            real_copy2 = shutil.copy2

            with patch("airr_igblast_paired.igblast.shutil.copy2", wraps=real_copy2) as copy2:
                _copy_file_if_needed(source, target)

            self.assertEqual(target.read_bytes(), b"new!")
            copy2.assert_called_once()
            temp_path = Path(copy2.call_args.args[1])
            self.assertEqual(temp_path.parent.resolve(), target.parent.resolve())
            self.assertNotEqual(temp_path, target)
            self.assertFalse(temp_path.exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_copy_file_reuses_same_content_without_copying(self) -> None:
        root = Path(f"test_igblast_copy_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            source = root / "source.bin"
            target = root / "target.bin"
            source.write_bytes(b"same content")
            target.write_bytes(b"same content")

            with patch("airr_igblast_paired.igblast.shutil.copy2") as copy2:
                _copy_file_if_needed(source, target)

            copy2.assert_not_called()
            self.assertEqual(target.read_bytes(), b"same content")
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_copy_file_failure_preserves_old_target_and_cleans_temp(self) -> None:
        root = Path(f"test_igblast_copy_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            source = root / "source.bin"
            target = root / "target.bin"
            source.write_bytes(b"new content")
            target.write_bytes(b"old content")

            def fail_copy(_: Path, temp_path: Path) -> None:
                Path(temp_path).write_bytes(b"partial")
                raise OSError("simulated copy failure")

            with patch("airr_igblast_paired.igblast.shutil.copy2", side_effect=fail_copy):
                with self.assertRaisesRegex(OSError, "simulated copy failure"):
                    _copy_file_if_needed(source, target)

            self.assertEqual(target.read_bytes(), b"old content")
            self.assertEqual({path.name for path in root.iterdir()}, {source.name, target.name})
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_copy_file_verification_failure_preserves_old_target_and_cleans_temp(self) -> None:
        root = Path(f"test_igblast_copy_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            source = root / "source.bin"
            target = root / "target.bin"
            source.write_bytes(b"new content")
            target.write_bytes(b"old content")

            def corrupt_copy(_: Path, temp_path: Path) -> None:
                Path(temp_path).write_bytes(b"bad content")

            with patch("airr_igblast_paired.igblast.shutil.copy2", side_effect=corrupt_copy):
                with self.assertRaisesRegex(OSError, "Copied file verification failed"):
                    _copy_file_if_needed(source, target)

            self.assertEqual(target.read_bytes(), b"old content")
            self.assertEqual({path.name for path in root.iterdir()}, {source.name, target.name})
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_default_threads_is_four(self) -> None:
        command = build_igblast_command(
            Path("queries.fasta"),
            Path("out.tsv"),
            IgBlastConfig(
                igblastn="igblastn",
                germline_db_v="human_gl_V",
                germline_db_j="human_gl_J",
            ),
        )

        self.assertEqual(command[command.index("-num_threads") + 1], "4")

    def test_build_igblast_command_uses_airr_outfmt_19(self) -> None:
        command = build_igblast_command(
            Path("queries.fasta"),
            Path("out.tsv"),
            IgBlastConfig(
                igblastn="igblastn",
                germline_db_v="human_gl_V",
                germline_db_d="human_gl_D",
                germline_db_j="human_gl_J",
                auxiliary_data="human_gl.aux",
                num_threads=4,
            ),
        )

        self.assertIn("-outfmt", command)
        self.assertEqual(command[command.index("-outfmt") + 1], "19")
        self.assertIn("-germline_db_D", command)
        self.assertIn("-auxiliary_data", command)
        self.assertEqual(command[command.index("-num_threads") + 1], "4")

    def test_subprocess_run_options_hide_windows_console(self) -> None:
        with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
            "airr_igblast_paired.igblast.subprocess.CREATE_NO_WINDOW",
            0x08000000,
            create=True,
        ):
            options = _subprocess_run_options()

        self.assertTrue(options["capture_output"])
        self.assertTrue(options["text"])
        self.assertFalse(options["check"])
        self.assertEqual(options["creationflags"], 0x08000000)

    def test_windows_runtime_context_prefers_igblast_internal_data(self) -> None:
        root = Path(f"test_igblast_runtime_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            exe = root / "tools" / "igblast-1.21.0" / "bin" / "igblastn.exe"
            exe.parent.mkdir(parents=True)
            exe.write_text("", encoding="utf-8")
            install_internal_data = exe.parent.parent / "internal_data"
            (install_internal_data / "human").mkdir(parents=True)
            (install_internal_data / "human" / "human_V.nsq").write_text("", encoding="utf-8")
            runtime_root = root / "runtime"

            refdata_root = root / "refdata" / "IgBlast_refdata_edit_imgt"
            refdata_internal_data = refdata_root / "internal_data"
            refdata_internal_data.mkdir(parents=True)
            db_dir = refdata_root / "db"
            db_dir.mkdir()
            v_db = db_dir / "IMGT_IGHV.imgt"
            j_db = db_dir / "IMGT_IGHJ.imgt"

            command = [
                str(exe),
                "-germline_db_V",
                str(v_db),
                "-germline_db_J",
                str(j_db),
            ]

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_short_path",
                side_effect=lambda path: str(path),
            ), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                cwd, env = _igblast_runtime_context(command)

            self.assertIsNotNone(cwd)
            assert cwd is not None
            self.assertEqual(cwd.parent, runtime_root / "internal_data_bundles")
            self.assertEqual(env["IGDATA"], ".")
            self.assertTrue((cwd / "human" / "human_V.nsq").exists())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_internal_data_bundles_are_immutable_per_install_and_content(self) -> None:
        root = Path(f"test_igblast_runtime_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            commands: list[list[str]] = []
            for install_name, content in (("install_a", "internal A"), ("install_b", "internal B")):
                exe = root / install_name / "bin" / "igblastn.exe"
                exe.parent.mkdir(parents=True)
                exe.write_text("", encoding="utf-8")
                internal_file = exe.parent.parent / "internal_data" / "human" / "human_V.nsq"
                internal_file.parent.mkdir(parents=True)
                internal_file.write_text(content, encoding="utf-8")
                commands.append([str(exe), "-germline_db_V", "v", "-germline_db_J", "j"])
            runtime_root = root / "runtime"

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_short_path",
                side_effect=lambda path: str(path),
            ), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                first_cwd, _ = _igblast_runtime_context(commands[0])
                second_cwd, _ = _igblast_runtime_context(commands[1])
                with patch("airr_igblast_paired.igblast.shutil.copy2") as copy2:
                    first_again, _ = _igblast_runtime_context(commands[0])

            self.assertIsNotNone(first_cwd)
            self.assertIsNotNone(second_cwd)
            assert first_cwd is not None and second_cwd is not None
            self.assertNotEqual(first_cwd, second_cwd)
            self.assertEqual(first_cwd, first_again)
            copy2.assert_not_called()
            self.assertEqual(
                (first_cwd / "human" / "human_V.nsq").read_text(encoding="utf-8"),
                "internal A",
            )
            self.assertEqual(
                (second_cwd / "human" / "human_V.nsq").read_text(encoding="utf-8"),
                "internal B",
            )
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_stages_germline_db_and_auxiliary_file(self) -> None:
        root = Path(f"test_igblast_stage_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            refdata_root = root / "refdata" / "IgBlast_refdata_edit_imgt"
            db_dir = refdata_root / "db"
            db_dir.mkdir(parents=True)
            prefix_name = "IMGT_IGHV.imgt"
            component_names = [
                f"{prefix_name}.ndb",
                f"{prefix_name}.not",
            ]
            (db_dir / component_names[0]).write_text("db index", encoding="utf-8")
            (db_dir / component_names[1]).write_text("db offsets", encoding="utf-8")
            (db_dir / f"{prefix_name}.njs").write_text(
                json.dumps(
                    {
                        "version": "1.2",
                        "dbtype": "Nucleotide",
                        "db-version": 5,
                        "files": component_names,
                    }
                ),
                encoding="utf-8",
            )
            optional_dir = refdata_root / "optional_file"
            optional_dir.mkdir()
            aux = optional_dir / "human_gl.aux"
            aux.write_text("aux", encoding="utf-8")
            runtime_root = root / "runtime"

            command = [
                "igblastn",
                "-germline_db_V",
                str(db_dir / prefix_name),
                "-auxiliary_data",
                str(aux),
            ]

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                staged = _stage_windows_igblast_resources(command)

            staged_v = Path(staged[staged.index("-germline_db_V") + 1])
            staged_aux = Path(staged[staged.index("-auxiliary_data") + 1])
            bundle_root = staged_v.parent.parent.parent
            self.assertEqual(bundle_root.parent, runtime_root / "resource_bundles")
            self.assertEqual(staged_v.parent.name, "V")
            self.assertEqual(staged_v.name, prefix_name)
            self.assertEqual(staged_aux, bundle_root / "optional_file" / "auxiliary_data.aux")
            staged_metadata = json.loads(
                Path(str(staged_v) + ".njs").read_text(encoding="utf-8")
            )
            self.assertEqual(staged_metadata["files"], component_names)
            for component_name in staged_metadata["files"]:
                self.assertTrue(
                    (staged_v.parent / component_name).is_file(),
                    component_name,
                )
            self.assertEqual(staged_aux.read_text(encoding="utf-8"), "aux")
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_resource_bundles_separate_same_basenames_by_source_and_content(self) -> None:
        root = Path(f"test_igblast_bundle_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            runtime_root = root / "runtime"

            def make_command(source_name: str, db_content: str, aux_content: str) -> list[str]:
                source_root = root / source_name
                db_dir = source_root / "db"
                db_dir.mkdir(parents=True)
                (db_dir / "shared.nsq").write_text(db_content, encoding="utf-8")
                optional_dir = source_root / "optional_file"
                optional_dir.mkdir()
                aux = optional_dir / "human_gl.aux"
                aux.write_text(aux_content, encoding="utf-8")
                return [
                    "igblastn",
                    "-germline_db_V",
                    str(db_dir / "shared"),
                    "-auxiliary_data",
                    str(aux),
                ]

            first_command = make_command("source_a", "database A", "aux A")
            second_command = make_command("source_b", "database B", "aux B")
            same_content_other_source = make_command("source_c", "database A", "aux A")

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                first = _stage_windows_igblast_resources(first_command)
                second = _stage_windows_igblast_resources(second_command)
                third = _stage_windows_igblast_resources(same_content_other_source)

            first_v = Path(first[first.index("-germline_db_V") + 1])
            second_v = Path(second[second.index("-germline_db_V") + 1])
            third_v = Path(third[third.index("-germline_db_V") + 1])
            first_bundle = first_v.parent.parent.parent
            second_bundle = second_v.parent.parent.parent
            third_bundle = third_v.parent.parent.parent
            self.assertEqual(len({first_bundle, second_bundle, third_bundle}), 3)
            self.assertEqual(Path(str(first_v) + ".nsq").read_text(encoding="utf-8"), "database A")
            self.assertEqual(Path(str(second_v) + ".nsq").read_text(encoding="utf-8"), "database B")
            self.assertEqual(Path(str(third_v) + ".nsq").read_text(encoding="utf-8"), "database A")
            self.assertEqual(
                (first_bundle / "optional_file" / "auxiliary_data.aux").read_text(encoding="utf-8"),
                "aux A",
            )
            self.assertEqual(
                (second_bundle / "optional_file" / "auxiliary_data.aux").read_text(encoding="utf-8"),
                "aux B",
            )
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_resource_bundle_reuses_same_immutable_paths(self) -> None:
        root = Path(f"test_igblast_bundle_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            db_dir = root / "source" / "db"
            db_dir.mkdir(parents=True)
            (db_dir / "shared.nsq").write_text("database", encoding="utf-8")
            optional_dir = root / "source" / "optional_file"
            optional_dir.mkdir()
            aux = optional_dir / "human_gl.aux"
            aux.write_text("aux", encoding="utf-8")
            command = [
                "igblastn",
                "-germline_db_V",
                str(db_dir / "shared"),
                "-auxiliary_data",
                str(aux),
            ]
            runtime_root = root / "runtime"

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                first = _stage_windows_igblast_resources(command)
                with patch("airr_igblast_paired.igblast.shutil.copy2") as copy2:
                    second = _stage_windows_igblast_resources(command)

            self.assertEqual(first, second)
            copy2.assert_not_called()
            self.assertEqual(len(list((runtime_root / "resource_bundles").iterdir())), 1)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_resource_source_mutation_does_not_publish_under_old_fingerprint(self) -> None:
        root = Path(f"test_igblast_bundle_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            db_dir = root / "source" / "db"
            db_dir.mkdir(parents=True)
            component = db_dir / "shared.nsq"
            component.write_text("old snapshot", encoding="utf-8")
            command = ["igblastn", "-germline_db_V", str(db_dir / "shared")]
            runtime_root = root / "runtime"
            real_copy2 = shutil.copy2

            def mutate_before_copy(source: Path, target: Path) -> None:
                Path(source).write_text("new snapshot", encoding="utf-8")
                real_copy2(source, target)

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ), patch(
                "airr_igblast_paired.igblast.shutil.copy2",
                side_effect=mutate_before_copy,
            ):
                with self.assertRaisesRegex(OSError, "source changed while staging"):
                    _stage_windows_igblast_resources(command)

            published_components = list(
                (runtime_root / "resource_bundles").glob("*/db/V/shared.nsq")
            )
            self.assertEqual(published_components, [])
            self.assertEqual(
                list((runtime_root / "resource_bundles").rglob("*.tmp")),
                [],
            )
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_windows_bundle_separates_v_d_j_with_the_same_prefix_basename(self) -> None:
        root = Path(f"test_igblast_bundle_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            command = ["igblastn"]
            for role in ("V", "D", "J"):
                db_dir = root / role / "db"
                db_dir.mkdir(parents=True)
                (db_dir / "shared.nsq").write_text(f"database {role}", encoding="utf-8")
                command.extend([f"-germline_db_{role}", str(db_dir / "shared")])
            runtime_root = root / "runtime"

            with patch("airr_igblast_paired.igblast.os.name", "nt"), patch(
                "airr_igblast_paired.igblast._windows_runtime_root",
                return_value=runtime_root,
            ):
                staged = _stage_windows_igblast_resources(command)

            staged_prefixes = {
                role: Path(staged[staged.index(f"-germline_db_{role}") + 1])
                for role in ("V", "D", "J")
            }
            self.assertEqual(len(set(staged_prefixes.values())), 3)
            self.assertEqual(
                len({path.parent.parent.parent for path in staged_prefixes.values()}),
                1,
            )
            for role, path in staged_prefixes.items():
                self.assertEqual(path.parent.name, role)
                self.assertEqual(path.name, "shared")
                self.assertEqual(
                    Path(str(path) + ".nsq").read_text(encoding="utf-8"),
                    f"database {role}",
                )
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_run_igblast_adds_disk_hint_for_iostream_error(self) -> None:
        with patch(
            "airr_igblast_paired.igblast._igblast_runtime_context",
            return_value=(None, {}),
        ), patch("airr_igblast_paired.igblast.subprocess.run") as fake_run:
            fake_run.return_value = SimpleNamespace(
                returncode=3,
                stdout="",
                stderr="WORKER: T2 BATCH # 256 EXCEPTION: ios_base::badbit set: iostream stream error",
            )

            with self.assertRaises(RuntimeError) as raised:
                run_igblast(
                    Path("queries.fasta"),
                    Path("out.airr.tsv"),
                    IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                )

        message = str(raised.exception)
        self.assertIn("Likely output/write failure", message)
        self.assertIn("available disk space", message)

    def test_run_igblast_batched_appends_one_header(self) -> None:
        root = Path(f"test_igblast_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            query = root / "queries.fasta"
            output = root / "out.airr.tsv"
            query.write_text(
                ">read1|R1\nAAA\n"
                ">read1|R2\nTTT\n"
                ">read2|R1\nCCC\n",
                encoding="utf-8",
            )

            def fake_run_igblast(
                batch_query: Path,
                batch_output: Path,
                _: IgBlastConfig,
                *,
                _resources_staged: bool = False,
                _runtime_context_override: tuple[Path | None, dict[str, str]] | None = None,
            ) -> list[str]:
                self.assertTrue(_resources_staged)
                self.assertIsNotNone(_runtime_context_override)
                names = [
                    line[1:].strip()
                    for line in batch_query.read_text(encoding="utf-8").splitlines()
                    if line.startswith(">")
                ]
                batch_output.write_text(
                    "sequence_id\tproductive\n"
                    + "".join(f"{name}\tT\n" for name in names),
                    encoding="utf-8",
                )
                return ["igblastn", "-query", str(batch_query), "-out", str(batch_output)]

            with patch(
                "airr_igblast_paired.igblast._stage_windows_igblast_config",
                wraps=_stage_windows_igblast_config,
            ) as stage_config, patch(
                "airr_igblast_paired.igblast._igblast_runtime_context",
                return_value=(None, {}),
            ) as runtime_context, patch(
                "airr_igblast_paired.igblast.run_igblast",
                side_effect=fake_run_igblast,
            ):
                command = run_igblast_batched(
                    query,
                    output,
                    IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                    batch_size=2,
                )

            lines = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines.count("sequence_id\tproductive"), 1)
            self.assertEqual(lines[1:], ["read1|R1\tT", "read1|R2\tT", "read2|R1\tT"])
            self.assertIn("# batches", command)
            self.assertIn("2", command)
            stage_config.assert_called_once()
            runtime_context.assert_called_once()
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
