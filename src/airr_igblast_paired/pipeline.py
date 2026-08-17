from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from contextlib import contextmanager
from collections.abc import Callable
import ctypes
import json
import os
import re
import shutil
import tempfile
import uuid

from . import __version__
from .igblast import IgBlastConfig, run_igblast, run_igblast_batched, validate_extra_igblast_args
from .pair_summary import PairSummaryStats, default_derived_tsv_paths, split_and_integrate_airr_tsv
from .prepare import (
    PrepareStats,
    ReadSelection,
    ReadTransform,
    ensure_distinct_paths,
    prepare_paired_fastq_to_fasta,
)
from .umi import UmiMode


CPM_MANIFEST_SCHEMA_VERSION = 2
CPM_COUNTING_SEMANTICS = "cpm_v3_read_pair_and_exact_raw_umi_per_clonotype_v1"
CPM_ANALYSIS_LABEL = "Read-pair annotation and exact UMI counting"

_ANALYSIS_SUFFIXES = ("umiSeq5", "umiCollapse1", "umiCollapse2", "umiNoCollapse", "noUmi")
_DYNAMIC_ANALYSIS_SUFFIX_RE = re.compile(r"_(?:umiSeq\d+|umiStrict[12])$")


@dataclass(frozen=True)
class PipelineResult:
    stats: PrepareStats
    command: list[str]
    query_fasta: Path | None
    output_tsv: Path
    r1_tsv: Path | None = None
    r2_tsv: Path | None = None
    integrated_tsv: Path | None = None
    counts_tsv: Path | None = None
    counts_xlsx: Path | None = None
    final_productive_counts_tsv: Path | None = None
    final_productive_counts_xlsx: Path | None = None
    umi_counts_tsv: Path | None = None
    umi_counts_xlsx: Path | None = None
    final_productive_umi_counts_tsv: Path | None = None
    final_productive_umi_counts_xlsx: Path | None = None
    pair_summary_stats: PairSummaryStats | None = None


@dataclass(frozen=True)
class NamedPipelineResult:
    label: str
    analysis_suffix: str
    result: PipelineResult


@dataclass(frozen=True)
class MultiPipelineResult:
    runs: tuple[NamedPipelineResult, ...]
    manifest_path: Path


@dataclass(frozen=True)
class _InputSnapshot:
    requested_path: Path
    resolved_path: Path
    size: int
    mtime_ns: int
    st_dev: int | None
    st_ino: int | None


def _strip_known_analysis_suffix(stem: str) -> str:
    for suffix in _ANALYSIS_SUFFIXES:
        marker = f"_{suffix}"
        if stem.endswith(marker):
            return stem[: -len(marker)]
    dynamic_match = _DYNAMIC_ANALYSIS_SUFFIX_RE.search(stem)
    if dynamic_match:
        return stem[: dynamic_match.start()]
    return stem


def path_with_analysis_suffix(path: str | Path, analysis_suffix: str) -> Path:
    path = Path(path)
    lower_name = path.name.lower()
    known_endings = (".airr.tsv", ".queries.fasta", ".queries.fa", ".fasta", ".fa", ".tsv")
    for ending in known_endings:
        if lower_name.endswith(ending):
            stem = path.name[: -len(ending)]
            base = _strip_known_analysis_suffix(stem)
            suffix = f"_{analysis_suffix}" if analysis_suffix else ""
            return path.with_name(f"{base}{suffix}{ending}")

    base = _strip_known_analysis_suffix(path.stem)
    suffix = f"_{analysis_suffix}" if analysis_suffix else ""
    return path.with_name(f"{base}{suffix}{path.suffix}")


def cpm_run_manifest_path(output_tsv: str | Path) -> Path:
    """Return the shared completion-manifest path for all aliases of a CPM run."""

    path = Path(output_tsv)
    lower_name = path.name.lower()
    for ending in (".airr.tsv", ".tsv"):
        if lower_name.endswith(ending):
            stem = path.name[: -len(ending)]
            break
    else:
        stem = path.stem
    base = _strip_known_analysis_suffix(stem)
    return path.with_name(f"{base}.run.json")


def _build_derived_outputs(
    output_tsv: Path,
) -> tuple[
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    Path,
    PairSummaryStats,
]:
    derived_paths, pair_stats = split_and_integrate_airr_tsv(output_tsv)
    return (
        derived_paths.r1_tsv,
        derived_paths.r2_tsv,
        derived_paths.integrated_tsv,
        derived_paths.counts_tsv,
        derived_paths.counts_xlsx,
        derived_paths.final_productive_counts_tsv,
        derived_paths.final_productive_counts_xlsx,
        derived_paths.umi_counts_tsv,
        derived_paths.umi_counts_xlsx,
        derived_paths.final_productive_umi_counts_tsv,
        derived_paths.final_productive_umi_counts_xlsx,
        pair_stats,
    )


def default_work_dir() -> Path:
    base = os.environ.get("LOCALAPPDATA")
    if base:
        return Path(base) / "PairedFastqIgblastAirrTsv" / "work"
    return Path(tempfile.gettempdir()) / "PairedFastqIgblastAirrTsv" / "work"


_RESULT_PATH_FIELDS = (
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
)


def planned_paired_output_paths(
    output_tsv: str | Path,
    query_fasta: str | Path | None = None,
) -> tuple[Path, ...]:
    output = Path(output_tsv)
    derived = default_derived_tsv_paths(output)
    paths = [
        output,
        derived.r1_tsv,
        derived.r2_tsv,
        derived.integrated_tsv,
        derived.counts_tsv,
        derived.counts_xlsx,
        derived.final_productive_counts_tsv,
        derived.final_productive_counts_xlsx,
        derived.umi_counts_tsv,
        derived.umi_counts_xlsx,
        derived.final_productive_umi_counts_tsv,
        derived.final_productive_umi_counts_xlsx,
    ]
    if query_fasta:
        paths.append(Path(query_fasta))
    return tuple(paths)


def _validate_run_paths(
    r1_path: str | Path,
    r2_path: str | Path,
    output_paths: tuple[Path, ...],
    *,
    igblast_config: IgBlastConfig | None = None,
) -> None:
    named_paths: dict[str, str | Path] = {
        "R1 FASTQ": r1_path,
        "R2 FASTQ": r2_path,
    }
    if igblast_config is not None:
        named_paths.update(_igblast_protected_paths(igblast_config))
    named_paths.update({f"output {index + 1}": path for index, path in enumerate(output_paths)})
    ensure_distinct_paths(**named_paths)


def _igblast_protected_paths(config: IgBlastConfig) -> dict[str, Path]:
    """Return IgBLAST executables/reference files that outputs must never replace."""

    protected: dict[str, Path] = {}
    executable_text = str(config.igblastn).strip()
    if executable_text:
        executable = Path(executable_text).expanduser()
        if not executable.is_file():
            resolved_executable = shutil.which(executable_text)
            executable = Path(resolved_executable) if resolved_executable else executable
        if executable.is_file():
            protected["IgBLAST executable"] = executable

    for label, prefix_text in (
        ("V DB", config.germline_db_v),
        ("D DB", config.germline_db_d),
        ("J DB", config.germline_db_j),
    ):
        if not prefix_text:
            continue
        prefix = Path(prefix_text).expanduser()
        protected[f"{label} prefix"] = prefix
        parent = prefix.parent
        if not parent.is_dir():
            continue
        marker = (prefix.name + ".").casefold()
        component_index = 0
        for candidate in parent.iterdir():
            if candidate.is_file() and candidate.name.casefold().startswith(marker):
                component_index += 1
                protected[f"{label} component {component_index}"] = candidate

    if config.auxiliary_data:
        protected["IgBLAST auxiliary data"] = Path(config.auxiliary_data).expanduser()
    return protected


def _validate_overwrite_policy(output_paths: tuple[Path, ...], *, overwrite: bool) -> None:
    invalid = [path for path in output_paths if path.exists() and not path.is_file()]
    if invalid:
        raise ValueError(
            "Output targets must be files, not directories or other filesystem objects:\n"
            + "\n".join(f"- {path}" for path in invalid)
        )
    if overwrite:
        return
    existing = [path for path in output_paths if path.exists()]
    if not existing:
        return
    preview = "\n".join(f"- {path}" for path in existing[:10])
    if len(existing) > 10:
        preview += f"\n- ... and {len(existing) - 10} more"
    raise FileExistsError(
        "Refusing to overwrite existing analysis outputs. "
        "Choose another output name/folder or explicitly enable overwrite:\n"
        f"{preview}"
    )


@contextmanager
def _output_locks(output_paths: tuple[Path, ...]):
    """Lock every final target so overlapping output sets cannot run together."""

    unique_targets: dict[str, Path] = {}
    for target in output_paths:
        resolved = target.expanduser().resolve(strict=False)
        unique_targets.setdefault(os.path.normcase(str(resolved)), resolved)

    lock_paths = [
        target.with_name(f".{target.name}.analysis.lock")
        for _, target in sorted(unique_targets.items())
    ]
    acquired: list[Path] = []
    try:
        for lock_path in lock_paths:
            lock_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError as exc:
                raise RuntimeError(
                    f"Another analysis may already be writing an overlapping output: {lock_path}. "
                    "If no analysis is running, remove the stale *.analysis.lock files "
                    "for this output set and retry."
                ) from exc
            acquired.append(lock_path)
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                handle.write(f"pid={os.getpid()}\n")
        yield
    finally:
        for lock_path in reversed(acquired):
            lock_path.unlink(missing_ok=True)


def _result_paths(result: PipelineResult) -> tuple[Path, ...]:
    paths: list[Path] = []
    for field in _RESULT_PATH_FIELDS:
        value = getattr(result, field)
        if value is not None:
            paths.append(Path(value))
    if result.query_fasta is not None:
        paths.append(Path(result.query_fasta))
    return tuple(paths)


def _result_with_final_paths(
    result: PipelineResult,
    *,
    output_tsv: Path,
    query_fasta: Path | None,
) -> PipelineResult:
    derived = default_derived_tsv_paths(output_tsv)
    return PipelineResult(
        stats=result.stats,
        command=result.command,
        query_fasta=query_fasta,
        output_tsv=output_tsv,
        r1_tsv=derived.r1_tsv,
        r2_tsv=derived.r2_tsv,
        integrated_tsv=derived.integrated_tsv,
        counts_tsv=derived.counts_tsv,
        counts_xlsx=derived.counts_xlsx,
        final_productive_counts_tsv=derived.final_productive_counts_tsv,
        final_productive_counts_xlsx=derived.final_productive_counts_xlsx,
        umi_counts_tsv=derived.umi_counts_tsv,
        umi_counts_xlsx=derived.umi_counts_xlsx,
        final_productive_umi_counts_tsv=derived.final_productive_umi_counts_tsv,
        final_productive_umi_counts_xlsx=derived.final_productive_umi_counts_xlsx,
        pair_summary_stats=result.pair_summary_stats,
    )


def _utc_timestamp(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, UTC).isoformat().replace("+00:00", "Z")


def _capture_input_snapshot(path: str | Path) -> _InputSnapshot:
    requested = Path(os.path.abspath(str(Path(path).expanduser())))
    resolved = requested.resolve(strict=True)
    stat = resolved.stat()
    if not resolved.is_file():
        raise ValueError(f"Input FASTQ is not a file: {resolved}")
    return _InputSnapshot(
        requested_path=requested,
        resolved_path=resolved,
        size=stat.st_size,
        mtime_ns=stat.st_mtime_ns,
        st_dev=getattr(stat, "st_dev", None),
        st_ino=getattr(stat, "st_ino", None),
    )


def _assert_input_snapshot_unchanged(label: str, snapshot: _InputSnapshot) -> None:
    try:
        resolved = snapshot.requested_path.resolve(strict=True)
        stat = resolved.stat()
    except OSError as exc:
        raise RuntimeError(f"{label} FASTQ changed or became unavailable during analysis") from exc

    identity = (
        os.path.normcase(str(resolved)),
        stat.st_size,
        stat.st_mtime_ns,
        getattr(stat, "st_dev", None),
        getattr(stat, "st_ino", None),
    )
    expected = (
        os.path.normcase(str(snapshot.resolved_path)),
        snapshot.size,
        snapshot.mtime_ns,
        snapshot.st_dev,
        snapshot.st_ino,
    )
    if identity != expected:
        raise RuntimeError(
            f"{label} FASTQ changed during analysis; no CPM output set was published"
        )


def _input_metadata(snapshot: _InputSnapshot) -> dict[str, str | int]:
    return {
        "path": str(snapshot.resolved_path),
        "size": snapshot.size,
        "mtime": _utc_timestamp(snapshot.mtime_ns / 1_000_000_000),
    }


def _manifest_output_entries(
    staged_result: PipelineResult,
    final_result: PipelineResult,
) -> list[dict[str, str | int]]:
    entries: list[dict[str, str | int]] = []
    for name in (*_RESULT_PATH_FIELDS, "query_fasta"):
        staged_value = getattr(staged_result, name)
        final_value = getattr(final_result, name)
        if staged_value is None and final_value is None:
            continue
        if staged_value is None or final_value is None:
            raise RuntimeError(f"staged/final output mismatch for manifest field: {name}")
        staged_path = Path(staged_value)
        if not staged_path.is_file():
            raise RuntimeError(f"Expected staged output is missing: {staged_path}")
        entries.append(
            {
                "name": name,
                "path": str(Path(final_value).expanduser().resolve(strict=False)),
                "size": staged_path.stat().st_size,
            }
        )
    return entries


def _write_cpm_run_manifest(
    path: Path,
    *,
    final_manifest_path: Path,
    r1_snapshot: _InputSnapshot,
    r2_snapshot: _InputSnapshot,
    read_selection: ReadSelection,
    r1_transform: ReadTransform,
    r2_transform: ReadTransform,
    query_name_template: str,
    min_length: int,
    max_n_rate: float,
    strict_ids: bool,
    umi_anchor_max_mismatches: int,
    igblast_batch_size: int | None,
    igblast_config: IgBlastConfig,
    staged_runs: list[NamedPipelineResult],
    final_runs: list[NamedPipelineResult],
) -> None:
    if len(staged_runs) != 1 or len(final_runs) != 1:
        raise RuntimeError("Ver3 manifest requires exactly one analysis result")

    modes = []
    for staged_run, final_run in zip(staged_runs, final_runs, strict=True):
        if staged_run.label != CPM_ANALYSIS_LABEL or final_run.label != CPM_ANALYSIS_LABEL:
            raise RuntimeError("run label mismatch while creating manifest")
        if staged_run.analysis_suffix or final_run.analysis_suffix:
            raise RuntimeError("run suffix mismatch while creating manifest")
        modes.append(
            {
                "label": CPM_ANALYSIS_LABEL,
                "analysis_suffix": "",
                "umi_counting": "exact_raw_umi_per_clonotype",
                "stats": {
                    "prepare": asdict(staged_run.result.stats),
                    "pair_summary": (
                        asdict(staged_run.result.pair_summary_stats)
                        if staged_run.result.pair_summary_stats is not None
                        else None
                    ),
                },
                "outputs": _manifest_output_entries(staged_run.result, final_run.result),
            }
        )

    payload = {
        "manifest_schema_version": CPM_MANIFEST_SCHEMA_VERSION,
        "software_version": __version__,
        "counting_semantics": CPM_COUNTING_SEMANTICS,
        "completed_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "manifest_path": str(final_manifest_path.expanduser().resolve(strict=False)),
        "inputs": {
            "r1": _input_metadata(r1_snapshot),
            "r2": _input_metadata(r2_snapshot),
        },
        "settings": {
            "run": {
                "read_selection": read_selection,
                "r1_transform": asdict(r1_transform),
                "r2_transform": asdict(r2_transform),
                "query_name_template": query_name_template,
                "strict_ids": strict_ids,
            },
            "qc": {
                "min_length": min_length,
                "max_n_rate": max_n_rate,
            },
            "umi": {
                "mode": "cpm-r2",
                "anchor_max_mismatches": umi_anchor_max_mismatches,
                "counting_unit": "exact_raw_umi_per_clonotype",
                "missing_umi": "retained_as_read_pair_support_not_umi_family",
            },
            "igblast": {
                "executable": igblast_config.igblastn,
                "germline_db_v": igblast_config.germline_db_v,
                "germline_db_d": igblast_config.germline_db_d,
                "germline_db_j": igblast_config.germline_db_j,
                "organism": igblast_config.organism,
                "domain_system": igblast_config.domain_system,
                "ig_seqtype": igblast_config.ig_seqtype,
                "auxiliary_data": igblast_config.auxiliary_data,
                "num_threads": igblast_config.num_threads,
                "extra_args": list(igblast_config.extra_args),
                "batch_size": igblast_batch_size,
            },
        },
        "modes": modes,
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _publish_files(
    staged_paths: tuple[Path, ...],
    final_paths: tuple[Path, ...],
    *,
    overwrite: bool,
    completion_marker: Path | None = None,
) -> None:
    if len(staged_paths) != len(final_paths):
        raise ValueError("staged/final output path count mismatch")
    if completion_marker is not None:
        completion_marker = Path(completion_marker)
        if completion_marker not in final_paths:
            raise ValueError("completion marker must be one of the final output paths")
        if final_paths[-1] != completion_marker:
            raise ValueError("completion marker must be the final published output")
    _validate_overwrite_policy(final_paths, overwrite=overwrite)

    transaction_id = uuid.uuid4().hex
    partials: list[tuple[Path, Path]] = []
    backups: list[tuple[Path, Path]] = []
    published: list[Path] = []
    try:
        for source, destination in zip(staged_paths, final_paths, strict=True):
            if not source.is_file():
                raise RuntimeError(f"Expected staged output is missing: {source}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            partial = destination.with_name(f".{destination.name}.{transaction_id}.partial")
            partials.append((partial, destination))
            shutil.copy2(source, partial)

        _validate_overwrite_policy(final_paths, overwrite=overwrite)
        if overwrite:
            backup_order = final_paths
            if completion_marker is not None:
                backup_order = (completion_marker,) + tuple(
                    destination for destination in final_paths if destination != completion_marker
                )
            for destination in backup_order:
                if not destination.exists():
                    continue
                backup = destination.with_name(f".{destination.name}.{transaction_id}.backup")
                backups.append((destination, backup))
                os.replace(destination, backup)

        for partial, destination in partials:
            published.append(destination)
            os.replace(partial, destination)
    except BaseException as publish_error:
        rollback_errors: list[str] = []
        published_attempts = set(published)

        # The completion marker is the transaction commit point and is always
        # published last.  If its atomic rename completed before an exception
        # (for example Ctrl+C delivered immediately after os.replace), every
        # data destination has already been published.  Rolling data back now
        # would risk leaving a new, locked marker beside old data.  Keep the
        # committed new set instead and only discard the obsolete backups.
        if completion_marker is not None and completion_marker in published_attempts:
            try:
                completion_marker.stat()
            except FileNotFoundError:
                pass
            except OSError as marker_state_error:
                raise RuntimeError(
                    "Output publishing was interrupted while the completion-marker "
                    "state could not be verified. Do not use or replace the output "
                    "set until it is inspected; hidden backups were retained. "
                    f"Marker check failed for {completion_marker}: {marker_state_error}"
                ) from publish_error
            else:
                for _, backup in backups:
                    try:
                        backup.unlink(missing_ok=True)
                    except OSError:
                        # The formal marker and data are the committed set. A
                        # leftover hidden backup is safe to inspect/remove later.
                        pass
                raise

        moved_backups: list[tuple[Path, Path]] = []
        for destination, backup in backups:
            if backup.exists():
                moved_backups.append((destination, backup))
            elif destination in published_attempts or not destination.exists():
                rollback_errors.append(f"recorded backup is missing for {destination}: {backup}")

        backed_up_destinations = {destination for destination, _ in moved_backups}
        marker_backup = next(
            (
                (destination, backup)
                for destination, backup in moved_backups
                if completion_marker is not None and destination == completion_marker
            ),
            None,
        )
        marker_publication_attempted = (
            completion_marker is not None and completion_marker in published_attempts
        )
        if completion_marker is not None and (
            marker_backup is not None or marker_publication_attempted
        ):
            try:
                completion_marker.unlink(missing_ok=True)
            except OSError as exc:
                rollback_errors.append(f"remove completion marker {completion_marker}: {exc}")

        for destination, backup in reversed(moved_backups):
            if completion_marker is not None and destination == completion_marker:
                continue
            try:
                os.replace(backup, destination)
            except OSError as exc:
                rollback_errors.append(f"restore {destination}: {exc}")
        for destination in reversed(published):
            if completion_marker is not None and destination == completion_marker:
                continue
            if destination in backed_up_destinations:
                continue
            try:
                destination.unlink(missing_ok=True)
            except OSError as exc:
                rollback_errors.append(f"remove new {destination}: {exc}")

        for destination, backup in backups:
            if (destination, backup) in moved_backups:
                continue
            if not destination.exists():
                rollback_errors.append(f"original output is missing after rollback: {destination}")

        if not rollback_errors and marker_backup is not None:
            marker_destination, marker_backup_path = marker_backup
            try:
                os.replace(marker_backup_path, marker_destination)
            except OSError as exc:
                rollback_errors.append(f"restore completion marker {marker_destination}: {exc}")

        if rollback_errors:
            if completion_marker is not None:
                try:
                    completion_marker.unlink(missing_ok=True)
                except OSError:
                    pass
            raise RuntimeError(
                "Output publishing failed and rollback was incomplete. "
                "Do not use the output set until it is inspected. "
                + "; ".join(rollback_errors)
            ) from publish_error
        raise
    else:
        for _, backup in backups:
            try:
                backup.unlink(missing_ok=True)
            except OSError:
                # The completed formal outputs are authoritative. A leftover
                # hidden backup can be cleaned manually without invalidating them.
                pass
    finally:
        for partial, _ in partials:
            try:
                partial.unlink(missing_ok=True)
            except OSError:
                pass


@contextmanager
def keep_windows_awake():
    if os.name != "nt":
        yield
        return

    es_continuous = 0x80000000
    es_system_required = 0x00000001
    es_awaymode_required = 0x00000040
    ctypes.windll.kernel32.SetThreadExecutionState(es_continuous | es_system_required | es_awaymode_required)
    try:
        yield
    finally:
        ctypes.windll.kernel32.SetThreadExecutionState(es_continuous)


def _run_igblast_maybe_batched(
    query_fasta: Path,
    output_tsv: Path,
    igblast_config: IgBlastConfig,
    igblast_batch_size: int | None,
    progress_callback: Callable[[str], None] | None = None,
) -> list[str]:
    if igblast_batch_size and igblast_batch_size > 0:
        return run_igblast_batched(
            query_fasta,
            output_tsv,
            igblast_config,
            batch_size=igblast_batch_size,
            progress_callback=progress_callback,
        )
    if progress_callback:
        progress_callback("Starting IgBLAST without batching...")
    return run_igblast(query_fasta, output_tsv, igblast_config)


def run_paired_igblast(
    *,
    r1_path: str | Path,
    r2_path: str | Path,
    output_tsv: str | Path,
    igblast_config: IgBlastConfig,
    query_fasta: str | Path | None = None,
    read_selection: ReadSelection = "both",
    r1_transform: ReadTransform = ReadTransform("forward"),
    r2_transform: ReadTransform = ReadTransform("reverse-complement"),
    query_name_template: str = "{read_id}|{read}",
    min_length: int = 0,
    max_n_rate: float = 1.0,
    strict_ids: bool = True,
    umi_mode: UmiMode = "none",
    umi_anchor_max_mismatches: int = 2,
    igblast_batch_size: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    work_dir: str | Path | None = None,
    overwrite: bool = False,
    _acquire_output_lock: bool = True,
) -> PipelineResult:
    if not str(r1_path).strip():
        raise ValueError("R1 FASTQ is required")
    if not str(r2_path).strip():
        raise ValueError("R2 FASTQ is required")
    if not str(output_tsv).strip():
        raise ValueError("Output TSV is required")
    if not igblast_config.germline_db_v.strip():
        raise ValueError("V DB prefix is required")
    if not igblast_config.germline_db_j.strip():
        raise ValueError("J DB prefix is required")
    if igblast_config.num_threads < 1:
        raise ValueError("IgBLAST thread count must be 1 or greater")
    if igblast_batch_size is not None and igblast_batch_size < 0:
        raise ValueError("IgBLAST batch size must be 0 or greater")
    validate_extra_igblast_args(igblast_config.extra_args)

    output_tsv = Path(output_tsv)
    final_query_fasta = Path(query_fasta) if query_fasta else None
    final_paths = planned_paired_output_paths(output_tsv, final_query_fasta)
    _validate_run_paths(
        r1_path,
        r2_path,
        final_paths,
        igblast_config=igblast_config,
    )
    _validate_overwrite_policy(final_paths, overwrite=overwrite)

    work_root = (Path(work_dir) if work_dir else default_work_dir()).expanduser().resolve(strict=False)
    work_root.mkdir(parents=True, exist_ok=True)
    scratch_dir = work_root / f"run.{uuid.uuid4().hex}"
    scratch_dir.mkdir()
    scratch_output = scratch_dir / "analysis.airr.tsv"
    scratch_query = scratch_dir / "queries.fasta"

    @contextmanager
    def maybe_lock():
        if _acquire_output_lock:
            with _output_locks(final_paths):
                yield
        else:
            yield

    try:
        with maybe_lock(), keep_windows_awake():
            _validate_overwrite_policy(final_paths, overwrite=overwrite)
            if progress_callback:
                progress_callback("Preparing IgBLAST query FASTA in scratch...")
            stats = prepare_paired_fastq_to_fasta(
                r1_path,
                r2_path,
                scratch_query,
                read_selection=read_selection,
                r1_transform=r1_transform,
                r2_transform=r2_transform,
                query_name_template=query_name_template,
                min_length=min_length,
                max_n_rate=max_n_rate,
                strict_ids=strict_ids,
                umi_mode=umi_mode,
                umi_anchor_max_mismatches=umi_anchor_max_mismatches,
            )
            command = _run_igblast_maybe_batched(
                scratch_query,
                scratch_output,
                igblast_config,
                igblast_batch_size,
                progress_callback,
            )
            if stats.records_written > 0 and (
                not scratch_output.is_file() or scratch_output.stat().st_size == 0
            ):
                raise RuntimeError("IgBLAST produced an empty AIRR TSV for a non-empty query FASTA")
            if progress_callback:
                progress_callback("Creating R1/R2, integrated, and count outputs in scratch...")
            (
                r1_tsv,
                r2_tsv,
                integrated_tsv,
                counts_tsv,
                counts_xlsx,
                final_productive_counts_tsv,
                final_productive_counts_xlsx,
                umi_counts_tsv,
                umi_counts_xlsx,
                final_productive_umi_counts_tsv,
                final_productive_umi_counts_xlsx,
                pair_stats,
            ) = _build_derived_outputs(scratch_output)
            observed_airr_rows = pair_stats.r1_rows + pair_stats.r2_rows
            if observed_airr_rows != stats.records_written:
                raise RuntimeError(
                    "IgBLAST AIRR row count does not match the query FASTA: "
                    f"queries={stats.records_written}, AIRR rows={observed_airr_rows}"
                )
            staged_result = PipelineResult(
                stats=stats,
                command=command,
                query_fasta=scratch_query if final_query_fasta else None,
                output_tsv=scratch_output,
                r1_tsv=r1_tsv,
                r2_tsv=r2_tsv,
                integrated_tsv=integrated_tsv,
                counts_tsv=counts_tsv,
                counts_xlsx=counts_xlsx,
                final_productive_counts_tsv=final_productive_counts_tsv,
                final_productive_counts_xlsx=final_productive_counts_xlsx,
                umi_counts_tsv=umi_counts_tsv,
                umi_counts_xlsx=umi_counts_xlsx,
                final_productive_umi_counts_tsv=final_productive_umi_counts_tsv,
                final_productive_umi_counts_xlsx=final_productive_umi_counts_xlsx,
                pair_summary_stats=pair_stats,
            )
            if progress_callback:
                progress_callback("Publishing completed outputs...")
            _publish_files(
                _result_paths(staged_result),
                final_paths,
                overwrite=overwrite,
            )
            return _result_with_final_paths(
                staged_result,
                output_tsv=output_tsv,
                query_fasta=final_query_fasta,
            )
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)


def planned_cpm_output_paths(
    output_tsv: str | Path,
    query_fasta: str | Path | None = None,
) -> tuple[Path, ...]:
    run_output = path_with_analysis_suffix(output_tsv, "")
    run_query = path_with_analysis_suffix(query_fasta, "") if query_fasta else None
    planned = list(planned_paired_output_paths(run_output, run_query))
    planned.append(cpm_run_manifest_path(output_tsv))
    return tuple(planned)


def run_cpm_umi_igblast_outputs(
    *,
    r1_path: str | Path,
    r2_path: str | Path,
    output_tsv: str | Path,
    igblast_config: IgBlastConfig,
    query_fasta: str | Path | None = None,
    read_selection: ReadSelection = "both",
    r1_transform: ReadTransform = ReadTransform("forward"),
    r2_transform: ReadTransform = ReadTransform("reverse-complement"),
    query_name_template: str = "{read_id}|{read}|UMI={umi}",
    min_length: int = 0,
    max_n_rate: float = 1.0,
    strict_ids: bool = True,
    umi_anchor_max_mismatches: int = 2,
    igblast_batch_size: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    work_dir: str | Path | None = None,
    overwrite: bool = False,
) -> MultiPipelineResult:
    if not str(igblast_config.germline_db_d or "").strip():
        raise ValueError("D DB prefix is required for CPM analysis")
    if not str(igblast_config.auxiliary_data or "").strip():
        raise ValueError("IgBLAST auxiliary data is required for CPM analysis")
    output_tsv = path_with_analysis_suffix(output_tsv, "")
    query_fasta = path_with_analysis_suffix(query_fasta, "") if query_fasta else None
    final_paths = planned_cpm_output_paths(output_tsv, query_fasta)
    final_manifest_path = cpm_run_manifest_path(output_tsv)
    if not final_paths or final_paths[-1] != final_manifest_path:
        raise RuntimeError("CPM completion manifest must be the final planned output")
    _validate_run_paths(
        r1_path,
        r2_path,
        final_paths,
        igblast_config=igblast_config,
    )
    _validate_overwrite_policy(final_paths, overwrite=overwrite)
    r1_snapshot = _capture_input_snapshot(r1_path)
    r2_snapshot = _capture_input_snapshot(r2_path)
    _validate_run_paths(
        r1_snapshot.resolved_path,
        r2_snapshot.resolved_path,
        final_paths,
        igblast_config=igblast_config,
    )

    output_tsv = Path(output_tsv)
    query_fasta = Path(query_fasta) if query_fasta else None
    work_root = (Path(work_dir) if work_dir else default_work_dir()).expanduser().resolve(strict=False)
    work_root.mkdir(parents=True, exist_ok=True)
    transaction_dir = work_root / f"cpm.{uuid.uuid4().hex}"
    transaction_dir.mkdir()

    staged_runs: list[NamedPipelineResult] = []
    final_runs: list[NamedPipelineResult] = []
    try:
        with _output_locks(final_paths):
            _validate_overwrite_policy(final_paths, overwrite=overwrite)
            if progress_callback:
                progress_callback(f"Starting {CPM_ANALYSIS_LABEL}...")
            _assert_input_snapshot_unchanged("R1", r1_snapshot)
            _assert_input_snapshot_unchanged("R2", r2_snapshot)
            staged_output = transaction_dir / "analysis.airr.tsv"
            staged_query = transaction_dir / "analysis.queries.fasta" if query_fasta else None
            staged_result = run_paired_igblast(
                r1_path=r1_path,
                r2_path=r2_path,
                output_tsv=staged_output,
                igblast_config=igblast_config,
                query_fasta=staged_query,
                read_selection=read_selection,
                r1_transform=r1_transform,
                r2_transform=r2_transform,
                query_name_template=query_name_template,
                min_length=min_length,
                max_n_rate=max_n_rate,
                strict_ids=strict_ids,
                umi_mode="cpm-r2",
                umi_anchor_max_mismatches=umi_anchor_max_mismatches,
                igblast_batch_size=igblast_batch_size,
                progress_callback=progress_callback,
                work_dir=work_root,
                overwrite=False,
                _acquire_output_lock=False,
            )
            _assert_input_snapshot_unchanged("R1", r1_snapshot)
            _assert_input_snapshot_unchanged("R2", r2_snapshot)
            staged_runs.append(NamedPipelineResult(CPM_ANALYSIS_LABEL, "", staged_result))
            final_result = _result_with_final_paths(
                staged_result,
                output_tsv=output_tsv,
                query_fasta=query_fasta,
            )
            final_runs.append(NamedPipelineResult(CPM_ANALYSIS_LABEL, "", final_result))

            staged_manifest_path = transaction_dir / "completion.run.json"
            _assert_input_snapshot_unchanged("R1", r1_snapshot)
            _assert_input_snapshot_unchanged("R2", r2_snapshot)
            _write_cpm_run_manifest(
                staged_manifest_path,
                final_manifest_path=final_manifest_path,
                r1_snapshot=r1_snapshot,
                r2_snapshot=r2_snapshot,
                read_selection=read_selection,
                r1_transform=r1_transform,
                r2_transform=r2_transform,
                query_name_template=query_name_template,
                min_length=min_length,
                max_n_rate=max_n_rate,
                strict_ids=strict_ids,
                umi_anchor_max_mismatches=umi_anchor_max_mismatches,
                igblast_batch_size=igblast_batch_size,
                igblast_config=igblast_config,
                staged_runs=staged_runs,
                final_runs=final_runs,
            )
            if progress_callback:
                progress_callback("Publishing completed CPM outputs...")
            _assert_input_snapshot_unchanged("R1", r1_snapshot)
            _assert_input_snapshot_unchanged("R2", r2_snapshot)
            staged_paths = tuple(
                path
                for named_result in staged_runs
                for path in _result_paths(named_result.result)
            ) + (staged_manifest_path,)
            _publish_files(
                staged_paths,
                final_paths,
                overwrite=overwrite,
                completion_marker=final_manifest_path,
            )
        return MultiPipelineResult(tuple(final_runs), final_manifest_path)
    finally:
        shutil.rmtree(transaction_dir, ignore_errors=True)
