from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import ctypes
import tempfile
from collections.abc import Callable, Iterator


_MANAGED_IGBLAST_FLAGS = {
    "-query",
    "-out",
    "-outfmt",
    "-organism",
    "-domain_system",
    "-ig_seqtype",
    "-num_threads",
    "-germline_db_v",
    "-germline_db_d",
    "-germline_db_j",
    "-auxiliary_data",
}


@dataclass(frozen=True)
class IgBlastConfig:
    germline_db_v: str
    germline_db_j: str
    germline_db_d: str | None = None
    igblastn: str = "igblastn"
    organism: str = "human"
    domain_system: str = "imgt"
    ig_seqtype: str = "Ig"
    auxiliary_data: str | None = None
    num_threads: int = 4
    extra_args: list[str] = field(default_factory=list)


def validate_extra_igblast_args(extra_args: list[str]) -> None:
    """Prevent user arguments from overriding pipeline-owned paths/settings."""

    for token in extra_args:
        flag = str(token).split("=", 1)[0].lower()
        if flag.startswith("--"):
            flag = "-" + flag[2:]
        if flag in _MANAGED_IGBLAST_FLAGS:
            raise ValueError(
                f"extra IgBLAST arguments may not override pipeline-managed flag {token!r}"
            )


def build_igblast_command(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
) -> list[str]:
    validate_extra_igblast_args(config.extra_args)
    command = [
        config.igblastn,
        "-query",
        str(query_fasta),
        "-out",
        str(output_tsv),
        "-outfmt",
        "19",
        "-organism",
        config.organism,
        "-domain_system",
        config.domain_system,
        "-ig_seqtype",
        config.ig_seqtype,
        "-num_threads",
        str(config.num_threads),
    ]

    if config.germline_db_v:
        command.extend(["-germline_db_V", config.germline_db_v])
    if config.germline_db_d:
        command.extend(["-germline_db_D", config.germline_db_d])
    if config.germline_db_j:
        command.extend(["-germline_db_J", config.germline_db_j])
    if config.auxiliary_data:
        command.extend(["-auxiliary_data", config.auxiliary_data])
    command.extend(config.extra_args)
    return command


def _windows_short_path(path: str | Path) -> str:
    text = str(path)
    if os.name != "nt" or not text:
        return text

    buffer = ctypes.create_unicode_buffer(32768)
    result = ctypes.windll.kernel32.GetShortPathNameW(str(text), buffer, len(buffer))
    if result:
        return buffer.value
    return text


def _db_prefix_to_windows_short_path(prefix: str) -> str:
    if os.name != "nt" or not prefix:
        return prefix

    path = Path(prefix)
    if path.exists():
        return _windows_short_path(path)

    db_suffixes = (
        ".ndb",
        ".nhr",
        ".nin",
        ".nog",
        ".nos",
        ".not",
        ".nsq",
        ".ntf",
        ".nto",
        ".phr",
        ".pin",
        ".pog",
        ".psd",
        ".psi",
        ".psq",
    )
    if any(Path(str(path) + suffix).exists() for suffix in db_suffixes):
        return str(Path(_windows_short_path(path.parent)) / path.name)
    return prefix


def _file_to_windows_short_path(path_text: str) -> str:
    if os.name != "nt" or not path_text:
        return path_text

    path = Path(path_text)
    if path.exists():
        return _windows_short_path(path)
    if path.parent.exists():
        return str(Path(_windows_short_path(path.parent)) / path.name)
    return path_text


def _normalize_command_for_windows(command: list[str]) -> list[str]:
    if os.name != "nt":
        return command

    normalized = list(command)
    path_flags = {
        "-query": "file",
        "-out": "file",
        "-germline_db_V": "db",
        "-germline_db_D": "db",
        "-germline_db_J": "db",
        "-auxiliary_data": "file",
    }
    for index, value in enumerate(normalized):
        if index == 0:
            resolved = shutil.which(value) or value
            normalized[index] = _windows_short_path(resolved) if Path(resolved).exists() else value
            continue
        previous = normalized[index - 1] if index > 0 else ""
        if previous in path_flags:
            if path_flags[previous] == "db":
                normalized[index] = _db_prefix_to_windows_short_path(value)
            else:
                normalized[index] = _file_to_windows_short_path(value)
    return normalized


def _command_value(command: list[str], flag: str) -> str | None:
    try:
        index = command.index(flag)
    except ValueError:
        return None
    if index + 1 >= len(command):
        return None
    return command[index + 1]


def _windows_runtime_root() -> Path:
    roots: list[Path] = []
    if os.environ.get("LOCALAPPDATA"):
        roots.append(Path(os.environ["LOCALAPPDATA"]) / "PairedFastqIgblastAirrTsv" / "igblast_runtime")
    if os.environ.get("PUBLIC"):
        roots.append(Path(os.environ["PUBLIC"]) / "PairedFastqIgblastAirrTsv" / "igblast_runtime")
    roots.append(Path(tempfile.gettempdir()) / "PairedFastqIgblastAirrTsv" / "igblast_runtime")

    for root in roots:
        try:
            root.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        return root

    fallback = Path(tempfile.gettempdir()) / "PairedFastqIgblastAirrTsv" / "igblast_runtime"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_file_if_needed(
    source: Path,
    target: Path,
    *,
    expected_size: int | None = None,
    expected_sha256: str | None = None,
) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    source_size = source.stat().st_size if expected_size is None else expected_size
    source_digest = expected_sha256
    if target.exists() and target.stat().st_size == source_size:
        if source_digest is None:
            source_digest = _sha256_file(source)
        if _sha256_file(target) == source_digest:
            return

    if source_digest is None:
        source_digest = _sha256_file(source)

    file_descriptor, temp_name = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=target.parent,
    )
    os.close(file_descriptor)
    temp_path = Path(temp_name)
    try:
        shutil.copy2(source, temp_path)
        if temp_path.stat().st_size != source_size or _sha256_file(temp_path) != source_digest:
            raise OSError(
                f"Copied file verification failed or source changed while staging: "
                f"{source} -> {target}"
            )
        os.replace(temp_path, target)
    finally:
        temp_path.unlink(missing_ok=True)


def _resolved_source_key(path: Path) -> str:
    return os.path.normcase(str(path.expanduser().resolve(strict=False)))


def _file_snapshot(
    source: Path,
    cache: dict[str, tuple[int, str]],
) -> tuple[int, str]:
    key = _resolved_source_key(source)
    prior = cache.get(key)
    if prior is not None:
        return prior

    before = source.stat()
    file_digest = _sha256_file(source)
    after = source.stat()
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or getattr(before, "st_ino", None) != getattr(after, "st_ino", None)
    ):
        raise OSError(f"IgBLAST resource changed while fingerprinting: {source}")
    snapshot = (after.st_size, file_digest)
    cache[key] = snapshot
    return snapshot


def _db_component_files(prefix: str) -> tuple[Path, tuple[Path, ...]]:
    path = Path(prefix).expanduser()
    parent = path.parent
    if not parent.is_dir():
        return path, ()

    marker = os.path.normcase(path.name + ".")
    files = tuple(
        sorted(
            (
                source
                for source in parent.iterdir()
                if source.is_file()
                and os.path.normcase(source.name).startswith(marker)
            ),
            key=lambda source: os.path.normcase(source.name),
        )
    )
    return path, files


def _resource_bundle_fingerprint(
    db_resources: tuple[tuple[str, Path, tuple[Path, ...]], ...],
    auxiliary: Path | None,
) -> tuple[str, dict[str, tuple[int, str]]]:
    """Hash source identity and content for one immutable IgBLAST resource set."""

    digest = hashlib.sha256()
    digest.update(b"cpm-igblast-resource-bundle-v1\0")
    snapshots: dict[str, tuple[int, str]] = {}
    for role, prefix, files in db_resources:
        digest.update(role.encode("ascii") + b"\0")
        digest.update(_resolved_source_key(prefix).encode("utf-8") + b"\0")
        for source in files:
            source_size, source_digest = _file_snapshot(source, snapshots)
            digest.update(_resolved_source_key(source).encode("utf-8") + b"\0")
            digest.update(str(source_size).encode("ascii") + b"\0")
            digest.update(source_digest.encode("ascii") + b"\0")

    digest.update(b"AUX\0")
    if auxiliary is not None:
        source_size, source_digest = _file_snapshot(auxiliary, snapshots)
        digest.update(_resolved_source_key(auxiliary).encode("utf-8") + b"\0")
        digest.update(str(source_size).encode("ascii") + b"\0")
        digest.update(source_digest.encode("ascii") + b"\0")
    return digest.hexdigest(), snapshots


def _stage_db_prefix_in_bundle(
    role: str,
    prefix: Path,
    files: tuple[Path, ...],
    bundle_root: Path,
    snapshots: dict[str, tuple[int, str]],
) -> str:
    if not files:
        return str(prefix)

    # BLAST v5 metadata (notably ``.njs``) embeds component filenames using
    # the original database prefix basename.  Keep that basename byte-for-byte
    # and use the role directory to prevent V/D/J collisions instead of
    # renaming the prefix itself.
    target_prefix = bundle_root / "db" / role / prefix.name
    for source in files:
        source_size, source_digest = snapshots[_resolved_source_key(source)]
        component_suffix = source.name[len(prefix.name) :]
        _copy_file_if_needed(
            source,
            target_prefix.with_name(target_prefix.name + component_suffix),
            expected_size=source_size,
            expected_sha256=source_digest,
        )
    return str(target_prefix)


def _stage_file_in_bundle(
    path: Path,
    subfolder: str,
    bundle_root: Path,
    snapshots: dict[str, tuple[int, str]],
) -> str:
    target = bundle_root / subfolder / "auxiliary_data.aux"
    source_size, source_digest = snapshots[_resolved_source_key(path)]
    _copy_file_if_needed(
        path,
        target,
        expected_size=source_size,
        expected_sha256=source_digest,
    )
    return str(target)


def _stage_windows_igblast_resources(command: list[str]) -> list[str]:
    """Stage one complete DB/aux set in a source/content-addressed bundle.

    Bundle paths are never shared by different source sets or different file
    contents. Concurrent analyses may therefore reuse the same immutable bytes,
    but cannot replace another analysis's V/D/J or auxiliary resources.
    """

    if os.name != "nt":
        return command

    staged = list(command)
    db_specs: list[tuple[str, str, Path, tuple[Path, ...]]] = []
    for role, flag in (
        ("V", "-germline_db_V"),
        ("D", "-germline_db_D"),
        ("J", "-germline_db_J"),
    ):
        value = _command_value(staged, flag)
        if value:
            prefix, files = _db_component_files(value)
            db_specs.append((role, flag, prefix, files))

    auxiliary: Path | None = None
    auxiliary_value = _command_value(staged, "-auxiliary_data")
    if auxiliary_value:
        candidate = Path(auxiliary_value).expanduser()
        if candidate.is_file():
            auxiliary = candidate

    fingerprint, snapshots = _resource_bundle_fingerprint(
        tuple((role, prefix, files) for role, _, prefix, files in db_specs),
        auxiliary,
    )
    bundle_root = _windows_runtime_root() / "resource_bundles" / fingerprint

    for role, flag, prefix, files in db_specs:
        staged[staged.index(flag) + 1] = _stage_db_prefix_in_bundle(
            role,
            prefix,
            files,
            bundle_root,
            snapshots,
        )
    if auxiliary is not None:
        staged[staged.index("-auxiliary_data") + 1] = _stage_file_in_bundle(
            auxiliary,
            "optional_file",
            bundle_root,
            snapshots,
        )
    return staged


def _stage_windows_igblast_config(config: IgBlastConfig) -> IgBlastConfig:
    """Stage a config once so all batches in one run share one immutable bundle."""

    if os.name != "nt":
        return config

    probe_command = _stage_windows_igblast_resources(
        build_igblast_command("bundle-query.fasta", "bundle-output.airr.tsv", config)
    )
    return replace(
        config,
        germline_db_v=_command_value(probe_command, "-germline_db_V") or config.germline_db_v,
        germline_db_d=_command_value(probe_command, "-germline_db_D") or config.germline_db_d,
        germline_db_j=_command_value(probe_command, "-germline_db_J") or config.germline_db_j,
        auxiliary_data=_command_value(probe_command, "-auxiliary_data") or config.auxiliary_data,
    )


def _stage_internal_data_bundle(source_dir: Path) -> Path:
    """Return an immutable, source/content-addressed internal_data tree."""

    source_dir = source_dir.expanduser().resolve(strict=True)
    files = tuple(
        sorted(
            (source for source in source_dir.rglob("*") if source.is_file()),
            key=lambda source: os.path.normcase(source.relative_to(source_dir).as_posix()),
        )
    )
    snapshots: dict[str, tuple[int, str]] = {}
    digest = hashlib.sha256()
    digest.update(b"cpm-igblast-internal-data-bundle-v1\0")
    digest.update(_resolved_source_key(source_dir).encode("utf-8") + b"\0")
    for source in files:
        relative = source.relative_to(source_dir)
        source_size, source_digest = _file_snapshot(source, snapshots)
        digest.update(relative.as_posix().encode("utf-8") + b"\0")
        digest.update(str(source_size).encode("ascii") + b"\0")
        digest.update(source_digest.encode("ascii") + b"\0")

    target_dir = _windows_runtime_root() / "internal_data_bundles" / digest.hexdigest()
    target_dir.mkdir(parents=True, exist_ok=True)
    for source in files:
        source_size, source_digest = snapshots[_resolved_source_key(source)]
        _copy_file_if_needed(
            source,
            target_dir / source.relative_to(source_dir),
            expected_size=source_size,
            expected_sha256=source_digest,
        )
    return target_dir


def _available_space_text(path: Path) -> str:
    try:
        target = path if path.is_dir() else path.parent
        usage = shutil.disk_usage(target)
    except OSError:
        return "available disk space could not be checked"
    return f"available disk space near output folder: {usage.free / (1024**3):.2f} GB"


def _looks_like_output_write_failure(message: str) -> bool:
    lowered = message.lower()
    return (
        "ios_base::badbit" in lowered
        or "iostream stream error" in lowered
        or "no space left" in lowered
        or "not enough space" in lowered
        or "disk" in lowered and "space" in lowered
    )


def _refdata_root_from_command(command: list[str]) -> Path | None:
    for flag in ("-germline_db_V", "-germline_db_D", "-germline_db_J"):
        value = _command_value(command, flag)
        if not value:
            continue
        prefix = Path(value)
        parent = prefix.parent
        if parent.name.lower() != "db":
            continue
        root = parent.parent
        if (root / "internal_data").exists():
            return root
    return None


def _igblast_runtime_context(command: list[str]) -> tuple[Path | None, dict[str, str]]:
    env = os.environ.copy()
    if os.name != "nt":
        return None, env

    resolved = shutil.which(command[0]) or command[0]
    exe = Path(resolved)
    install_root = exe.parent.parent if exe.exists() else None
    internal_data = install_root / "internal_data" if install_root else None
    if internal_data and internal_data.exists():
        staged_internal_data = _stage_internal_data_bundle(internal_data)
        # IgBLAST looks up annotation files relative to IGDATA.  On Windows,
        # using an ASCII-friendly scratch cwd plus IGDATA='.' avoids BLAST
        # database path issues when the portable folder contains Japanese text.
        env["IGDATA"] = "."
        return Path(_windows_short_path(staged_internal_data)), env

    refdata_root = _refdata_root_from_command(command)
    if refdata_root is not None:
        env["IGDATA"] = _windows_short_path(refdata_root)
        return None, env
    return None, env


def _subprocess_run_options() -> dict[str, object]:
    options: dict[str, object] = {"capture_output": True, "text": True, "check": False}
    if os.name == "nt":
        options["creationflags"] = subprocess.CREATE_NO_WINDOW
    return options


def run_igblast(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
    *,
    _resources_staged: bool = False,
    _runtime_context_override: tuple[Path | None, dict[str, str]] | None = None,
) -> list[str]:
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = build_igblast_command(query_fasta, output_tsv, config)
    if not _resources_staged:
        command = _stage_windows_igblast_resources(command)
    command = _normalize_command_for_windows(command)
    if _runtime_context_override is None:
        cwd, env = _igblast_runtime_context(command)
    else:
        cwd, env = _runtime_context_override
    result = subprocess.run(command, cwd=cwd, env=env, **_subprocess_run_options())
    if result.returncode != 0:
        command_text = " ".join(command)
        message = result.stderr.strip() or result.stdout.strip() or "IgBLAST failed without output"
        if _looks_like_output_write_failure(message):
            message = (
                f"{message}\n\n"
                "Likely output/write failure. Free disk space and rerun the analysis. "
                f"{_available_space_text(output_tsv)}."
            )
        raise RuntimeError(f"IgBLAST failed with exit code {result.returncode}\n{command_text}\n{message}")
    return command


def _read_fasta_records(path: str | Path) -> Iterator[list[str]]:
    record: list[str] = []
    with Path(path).open("rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith(">") and record:
                yield record
                record = []
            record.append(line)
    if record:
        yield record


def _write_fasta_batch(path: Path, records: list[list[str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.writelines(record)


def _append_airr_tsv_batch(final_tsv: Path, batch_tsv: Path, *, wrote_header: bool) -> bool:
    with batch_tsv.open("rt", encoding="utf-8", newline="") as source, final_tsv.open(
        "at" if wrote_header else "wt",
        encoding="utf-8",
        newline="",
    ) as target:
        for line_number, line in enumerate(source):
            if wrote_header and line_number == 0 and line.startswith("sequence_id\t"):
                continue
            target.write(line)
            if line_number == 0 and line.startswith("sequence_id\t"):
                wrote_header = True
    return wrote_header


def run_igblast_batched(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
    *,
    batch_size: int,
    progress_callback: Callable[[str], None] | None = None,
) -> list[str]:
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than 0")

    query_fasta = Path(query_fasta)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.unlink(missing_ok=True)
    staged_config = _stage_windows_igblast_config(config)
    context_command = _normalize_command_for_windows(
        build_igblast_command(query_fasta, output_tsv, staged_config)
    )
    staged_runtime_context = _igblast_runtime_context(context_command)

    commands: list[list[str]] = []
    records: list[list[str]] = []
    batch_index = 0
    wrote_header = False

    def flush_batch() -> None:
        nonlocal batch_index, wrote_header, records
        if not records:
            return
        batch_index += 1
        batch_query = output_tsv.parent / f"{output_tsv.stem}.batch{batch_index:04d}.queries.fasta"
        batch_output = output_tsv.parent / f"{output_tsv.stem}.batch{batch_index:04d}.airr.tsv"
        query_count = len(records)
        try:
            if progress_callback:
                progress_callback(f"Starting IgBLAST batch {batch_index} ({query_count} queries)...")
            _write_fasta_batch(batch_query, records)
            commands.append(
                run_igblast(
                    batch_query,
                    batch_output,
                    staged_config,
                    _resources_staged=True,
                    _runtime_context_override=staged_runtime_context,
                )
            )
            wrote_header = _append_airr_tsv_batch(output_tsv, batch_output, wrote_header=wrote_header)
            if progress_callback:
                progress_callback(f"Finished IgBLAST batch {batch_index}.")
        finally:
            batch_query.unlink(missing_ok=True)
            batch_output.unlink(missing_ok=True)
            records = []

    for record in _read_fasta_records(query_fasta):
        records.append(record)
        if len(records) >= batch_size:
            flush_batch()
    flush_batch()

    if not commands:
        output_tsv.write_text("", encoding="utf-8")
        return []

    first_command = list(commands[0])
    first_command.extend(["# batches", str(len(commands)), "# batch_size", str(batch_size)])
    return first_command
