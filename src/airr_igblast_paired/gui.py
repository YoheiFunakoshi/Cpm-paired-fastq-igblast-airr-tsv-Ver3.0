from __future__ import annotations

from pathlib import Path
import queue
import shutil
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from .igblast import IgBlastConfig
from .naming import (
    PROJECT_FOLDER_NAME,
    default_data_folder,
    default_output_tsv_path,
    default_query_fasta_path,
    default_results_folder,
)
from .pipeline import (
    _validate_run_paths,
    default_work_dir,
    planned_cpm_output_paths,
    run_cpm_umi_igblast_outputs,
)
from .prepare import ReadTransform


APP_TITLE = PROJECT_FOLDER_NAME


_BLAST_DB_COMPONENT_SUFFIXES = {
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
}

_INTEGER_FORM_FIELDS = (
    ("num_threads", "Threads"),
    ("igblast_batch_size", "IgBLAST batch size"),
    ("trim_left_r1", "Trim left R1"),
    ("trim_right_r1", "Trim right R1"),
    ("trim_left_r2", "Trim left R2"),
    ("trim_right_r2", "Trim right R2"),
    ("min_length", "Min length"),
    ("umi_anchor_max_mismatches", "UMI anchor mismatches"),
)


def _has_blast_db_component(prefix: str | Path) -> bool:
    prefix_path = Path(prefix).expanduser()
    parent = prefix_path.parent
    if not parent.is_dir():
        return False
    return any(
        candidate.is_file() and candidate.suffix.lower() in _BLAST_DB_COMPONENT_SUFFIXES
        for candidate in parent.glob(prefix_path.name + ".*")
    )


def _has_blast_db(prefix: Path) -> bool:
    # Keep the legacy automatic refdata selection criterion unchanged.
    return Path(str(prefix) + ".nsq").exists()


def _normalize_run_values(raw_values: dict[str, object]) -> dict[str, object]:
    """Convert editable GUI fields to ordinary Python values on the Tk thread."""

    values = dict(raw_values)
    for key, label in _INTEGER_FORM_FIELDS:
        try:
            values[key] = int(values[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{label} には整数を入力してください") from exc
    try:
        values["max_n_rate"] = float(values["max_n_rate"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Max N rate には数値を入力してください") from exc
    values["strict_ids"] = bool(values.get("strict_ids", False))
    nonnegative_fields = (
        ("igblast_batch_size", "IgBLAST batch size"),
        ("trim_left_r1", "Trim left R1"),
        ("trim_right_r1", "Trim right R1"),
        ("trim_left_r2", "Trim left R2"),
        ("trim_right_r2", "Trim right R2"),
        ("min_length", "Min length"),
        ("umi_anchor_max_mismatches", "UMI anchor mismatches"),
    )
    for key, label in nonnegative_fields:
        if int(values[key]) < 0:
            raise ValueError(f"{label} は0以上にしてください")
    if int(values["num_threads"]) < 1:
        raise ValueError("Threads は1以上にしてください")
    if not 0 <= float(values["max_n_rate"]) <= 1:
        raise ValueError("Max N rate は0以上1以下にしてください")
    return values


def _preflight_run_errors(
    values: dict[str, object],
    planned_outputs: tuple[Path, ...],
) -> list[str]:
    """Return user-facing validation errors without invoking Tk or starting analysis."""

    errors: list[str] = []
    r1_path = Path(str(values.get("r1", ""))).expanduser()
    r2_path = Path(str(values.get("r2", ""))).expanduser()

    if not r1_path.is_file():
        errors.append(f"R1 FASTQ がファイルとして見つかりません: {r1_path}")
    if not r2_path.is_file():
        errors.append(f"R2 FASTQ がファイルとして見つかりません: {r2_path}")
    try:
        _validate_run_paths(
            r1_path,
            r2_path,
            planned_outputs,
            igblast_config=IgBlastConfig(
                igblastn=str(values.get("igblastn", "")),
                germline_db_v=str(values.get("germline_db_v", "")),
                germline_db_d=str(values.get("germline_db_d", "")) or None,
                germline_db_j=str(values.get("germline_db_j", "")),
                auxiliary_data=str(values.get("auxiliary_data", "")) or None,
            ),
        )
    except ValueError as exc:
        errors.append(f"R1/R2または出力先のパスが重複しています: {exc}")

    igblastn_value = str(values.get("igblastn", "")).strip()
    igblastn_path = Path(igblastn_value).expanduser() if igblastn_value else None
    resolved_igblastn = shutil.which(igblastn_value) if igblastn_value else None
    explicit_igblastn_exists = igblastn_path is not None and igblastn_path.is_file()
    path_igblastn_exists = resolved_igblastn is not None and Path(resolved_igblastn).is_file()
    if not (explicit_igblastn_exists or path_igblastn_exists):
        errors.append(f"igblastn 実行ファイルが見つかりません（指定パスまたはPATH）: {igblastn_value}")

    for key, label in (
        ("germline_db_v", "V DB prefix"),
        ("germline_db_d", "D DB prefix"),
        ("germline_db_j", "J DB prefix"),
    ):
        prefix = str(values.get(key, "")).strip()
        if not prefix or not _has_blast_db_component(prefix):
            errors.append(f"{label} のBLAST DB構成ファイル（.nsq等）が見つかりません: {prefix}")

    auxiliary_value = str(values.get("auxiliary_data", "")).strip()
    if not auxiliary_value or not Path(auxiliary_value).expanduser().is_file():
        errors.append(f"Aux file がファイルとして見つかりません: {auxiliary_value}")

    return errors


def _find_preferred_refdata_root() -> Path | None:
    data_folder = default_data_folder()
    desktop = Path.home() / "Desktop"
    legacy_data_folder = desktop / "CPM Paired Fastq IgBLAST AIRR tsv"
    local_refdata = data_folder / "refdata" / "IgBlast_refdata_edit_imgt"
    legacy_local_refdata = legacy_data_folder / "refdata" / "IgBlast_refdata_edit_imgt"
    desktop_refdata = desktop / "IgBlast_refdata_edit_imgt"

    candidates = [
        local_refdata,
        legacy_local_refdata,
        desktop_refdata,
    ]
    for root in candidates:
        if (
            _has_blast_db(root / "db" / "IMGT_IGHV.imgt")
            and _has_blast_db(root / "db" / "IMGT_IGHD.imgt")
            and _has_blast_db(root / "db" / "IMGT_IGHJ.imgt")
            and (root / "optional_file" / "human_gl.aux").exists()
        ):
            return root
    return None


def _find_preferred_igblast_root() -> Path | None:
    data_folder = default_data_folder()
    candidates = [
        data_folder / "tools" / "igblast-1.21.0",
        data_folder / "igblast-1.21.0",
        Path("C:/Program Files/NCBI/igblast-1.21.0"),
    ]
    for root in candidates:
        if (root / "bin" / "igblastn.exe").exists():
            return root
    return None


def _find_preferred_igblastn() -> str:
    igblast_root = _find_preferred_igblast_root()
    if igblast_root is not None:
        return str(igblast_root / "bin" / "igblastn.exe")
    return shutil.which("igblastn") or "igblastn"


class App(ttk.Frame):
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master, padding=12)
        self.master = master
        self.messages: queue.Queue[tuple[str, str]] = queue.Queue()
        self.vars: dict[str, tk.Variable] = {}
        self.data_folder = default_data_folder()
        self.results_folder = default_results_folder(self.data_folder)
        self._auto_output_path = ""
        self._auto_query_fasta_path = ""
        self._output_overridden = False
        self._query_fasta_overridden = False
        self._setting_auto_path = False
        self.running = False
        self._build_variables()
        self._build_layout()
        self._attach_path_traces()
        self._update_default_result_paths()
        self.master.protocol("WM_DELETE_WINDOW", self._on_close)
        self._poll_messages()

    def _build_variables(self) -> None:
        data_folder = self.data_folder
        imgt_vdj = Path.home() / "Desktop" / "IgWork" / "IMGT_VDJ"
        refdata_root = _find_preferred_refdata_root()
        igblast_root = _find_preferred_igblast_root()
        igblastn = _find_preferred_igblastn()
        aux_file = (
            refdata_root / "optional_file" / "human_gl.aux"
            if refdata_root
            else (igblast_root / "optional_file" / "human_gl.aux" if igblast_root else None)
        )

        if refdata_root:
            germline_db_v = refdata_root / "db" / "IMGT_IGHV.imgt"
            germline_db_d = refdata_root / "db" / "IMGT_IGHD.imgt"
            germline_db_j = refdata_root / "db" / "IMGT_IGHJ.imgt"
        else:
            germline_db_v = imgt_vdj / "human_IGHV_IMGT"
            germline_db_d = imgt_vdj / "human_IGHD_IMGT"
            germline_db_j = imgt_vdj / "human_IGHJ_IMGT"

        defaults: dict[str, str | int | float | bool] = {
            "r1": "",
            "r2": "",
            "out": "",
            "query_fasta": "",
            "igblastn": igblastn,
            "germline_db_v": str(germline_db_v) if _has_blast_db(germline_db_v) else "",
            "germline_db_d": str(germline_db_d) if _has_blast_db(germline_db_d) else "",
            "germline_db_j": str(germline_db_j) if _has_blast_db(germline_db_j) else "",
            "auxiliary_data": str(aux_file) if aux_file and aux_file.exists() else "",
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
        for key, value in defaults.items():
            if isinstance(value, bool):
                self.vars[key] = tk.BooleanVar(value=value)
            elif isinstance(value, int):
                self.vars[key] = tk.IntVar(value=value)
            elif isinstance(value, float):
                self.vars[key] = tk.DoubleVar(value=value)
            else:
                self.vars[key] = tk.StringVar(value=value)

    def _build_layout(self) -> None:
        self.grid(sticky="nsew")
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        row = 0
        row = self._add_file_row(row, "R1 FASTQ", "r1", "open")
        row = self._add_file_row(row, "R2 FASTQ", "r2", "open")
        row = self._add_file_row(row, "Output TSV", "out", "save")
        row = self._add_file_row(row, "Keep query FASTA", "query_fasta", "save")

        separator = ttk.Separator(self)
        separator.grid(row=row, column=0, columnspan=3, sticky="ew", pady=8)
        row += 1

        row = self._add_file_row(row, "igblastn", "igblastn", "open")
        row = self._add_db_row(row, "V DB prefix", "germline_db_v")
        row = self._add_db_row(row, "D DB prefix", "germline_db_d")
        row = self._add_db_row(row, "J DB prefix", "germline_db_j")
        row = self._add_file_row(row, "Aux file", "auxiliary_data", "open")

        row = self._add_entry_row(row, "Organism", "organism")
        row = self._add_entry_row(row, "Domain system", "domain_system")
        row = self._add_entry_row(row, "Seq type", "ig_seqtype")
        row = self._add_spin_row(row, "Threads", "num_threads", 1, 128)
        row = self._add_spin_row(row, "IgBLAST batch size", "igblast_batch_size", 0, 1000000)

        separator = ttk.Separator(self)
        separator.grid(row=row, column=0, columnspan=3, sticky="ew", pady=8)
        row += 1

        row = self._add_combo_row(row, "Reads", "read_selection", ("both", "r1", "r2"))
        row = self._add_combo_row(row, "R1 orientation", "r1_orientation", ("forward", "reverse-complement"))
        row = self._add_combo_row(row, "R2 orientation", "r2_orientation", ("forward", "reverse-complement"))
        row = self._add_spin_row(row, "Trim left R1", "trim_left_r1", 0, 10000)
        row = self._add_spin_row(row, "Trim right R1", "trim_right_r1", 0, 10000)
        row = self._add_spin_row(row, "Trim left R2", "trim_left_r2", 0, 10000)
        row = self._add_spin_row(row, "Trim right R2", "trim_right_r2", 0, 10000)
        row = self._add_spin_row(row, "Min length", "min_length", 0, 100000)
        row = self._add_entry_row(row, "Max N rate", "max_n_rate")
        row = self._add_entry_row(row, "Query name", "query_name_template")
        row = self._add_spin_row(row, "UMI anchor mismatches", "umi_anchor_max_mismatches", 0, 20)

        strict_ids = ttk.Checkbutton(self, text="Require matching R1/R2 IDs", variable=self.vars["strict_ids"])
        strict_ids.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        rule_label = ttk.Label(
            self,
            text="Final call rule: v_call uses R2 first; CDR3/junction_aa uses canonical side first, then R1",
        )
        rule_label.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        workflow_label = ttk.Label(
            self,
            text=(
                "Ver3: retained R1/R2 pairs are annotated once; RG-compatible R1/R2 and "
                "integrated read-pair outputs are created first"
            ),
            wraplength=760,
        )
        workflow_label.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        umi_label = ttk.Label(
            self,
            text=(
                "Exact UMI count: identical raw UMIs are one family only within each "
                "V/J/canonical junction-AA clonotype; UMI-missing pairs are retained"
            ),
            wraplength=760,
        )
        umi_label.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        productive_label = ttk.Label(
            self,
            text="Productive output: standard counts keep all included CDR3; final_productive_counts keeps final_productive=True only",
        )
        productive_label.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        self.run_button = ttk.Button(self, text="Run", command=self._start_run)
        self.run_button.grid(row=row, column=1, sticky="w", pady=8)
        row += 1

        self.log = scrolledtext.ScrolledText(self, height=10, width=80)
        self.log.grid(row=row, column=0, columnspan=3, sticky="nsew")
        self.rowconfigure(row, weight=1)

    def _add_file_row(self, row: int, label: str, key: str, mode: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, sticky="ew", pady=3)
        ttk.Button(self, text="Browse", command=lambda: self._browse(key, mode)).grid(row=row, column=2, padx=(8, 0))
        return row + 1

    def _add_db_row(self, row: int, label: str, key: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, sticky="ew", pady=3)
        ttk.Button(self, text="Browse", command=lambda: self._browse_db_prefix(key)).grid(
            row=row,
            column=2,
            padx=(8, 0),
        )
        return row + 1

    def _add_entry_row(self, row: int, label: str, key: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
        return row + 1

    def _add_spin_row(self, row: int, label: str, key: str, start: int, end: int) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Spinbox(self, from_=start, to=end, textvariable=self.vars[key], width=12).grid(
            row=row,
            column=1,
            sticky="w",
            pady=3,
        )
        return row + 1

    def _add_combo_row(self, row: int, label: str, key: str, values: tuple[str, ...]) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        combo = ttk.Combobox(self, textvariable=self.vars[key], values=values, state="readonly")
        combo.grid(row=row, column=1, sticky="w", pady=3)
        return row + 1

    def _browse(self, key: str, mode: str) -> None:
        if mode == "save":
            path = filedialog.asksaveasfilename(**self._save_dialog_options(key))
        else:
            path = filedialog.askopenfilename(**self._open_dialog_options(key))
        if path:
            self.vars[key].set(path)
            if key in {"r1", "r2"}:
                self._update_default_result_paths()

    def _open_dialog_options(self, key: str) -> dict[str, object]:
        options: dict[str, object] = {}
        if key in {"r1", "r2"} and self.data_folder.exists():
            options["initialdir"] = str(self.data_folder)
            options["filetypes"] = [
                ("FASTQ files", "*.fastq *.fastq.gz *.fq *.fq.gz"),
                ("All files", "*.*"),
            ]
        return options

    def _save_dialog_options(self, key: str) -> dict[str, object]:
        options: dict[str, object] = {}
        if key == "out":
            suggested = self._suggested_output_path()
            options["defaultextension"] = ".tsv"
            options["filetypes"] = [("TSV files", "*.tsv"), ("All files", "*.*")]
        elif key == "query_fasta":
            suggested = self._suggested_query_fasta_path()
            options["defaultextension"] = ".fasta"
            options["filetypes"] = [("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        else:
            suggested = None

        initialdir = suggested.parent if suggested else self.results_folder
        try:
            initialdir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        options["initialdir"] = str(initialdir)
        if suggested:
            options["initialfile"] = suggested.name
        return options

    def _attach_path_traces(self) -> None:
        self.vars["r1"].trace_add("write", lambda *_: self._update_default_result_paths())
        self.vars["r2"].trace_add("write", lambda *_: self._update_default_result_paths())
        self.vars["out"].trace_add("write", lambda *_: self._mark_output_overridden())
        self.vars["query_fasta"].trace_add("write", lambda *_: self._mark_query_fasta_overridden())

    def _mark_output_overridden(self) -> None:
        if self._setting_auto_path:
            return
        value = str(self.vars["out"].get()).strip()
        self._output_overridden = bool(value and value != self._auto_output_path)

    def _mark_query_fasta_overridden(self) -> None:
        if self._setting_auto_path:
            return
        value = str(self.vars["query_fasta"].get()).strip()
        self._query_fasta_overridden = bool(value and value != self._auto_query_fasta_path)

    def _suggested_output_path(self) -> Path | None:
        r1 = str(self.vars["r1"].get()).strip()
        r2 = str(self.vars["r2"].get()).strip()
        if not r1:
            return None
        return default_output_tsv_path(r1, r2 or None, self.data_folder)

    def _suggested_query_fasta_path(self) -> Path | None:
        r1 = str(self.vars["r1"].get()).strip()
        r2 = str(self.vars["r2"].get()).strip()
        if not r1:
            return None
        return default_query_fasta_path(r1, r2 or None, self.data_folder)

    def _update_default_result_paths(self) -> None:
        suggested_output = self._suggested_output_path()
        suggested_query = self._suggested_query_fasta_path()
        if not suggested_output:
            return

        current_output = str(self.vars["out"].get()).strip()
        current_query = str(self.vars["query_fasta"].get()).strip()
        if (
            self._output_overridden
            and current_output != self._auto_output_path
            and Path(current_output).name.lower() != "result.airr.tsv"
        ):
            update_output = False
        else:
            update_output = True

        if self._query_fasta_overridden and current_query != self._auto_query_fasta_path:
            update_query = False
        else:
            update_query = suggested_query is not None

        self._setting_auto_path = True
        try:
            if update_output:
                self._auto_output_path = str(suggested_output)
                self.vars["out"].set(self._auto_output_path)
                self._output_overridden = False
            if update_query and suggested_query is not None:
                self._auto_query_fasta_path = str(suggested_query)
                self.vars["query_fasta"].set(self._auto_query_fasta_path)
                self._query_fasta_overridden = False
        finally:
            self._setting_auto_path = False

    def _browse_db_prefix(self, key: str) -> None:
        path = filedialog.askopenfilename()
        if not path:
            return
        selected = Path(path)
        if selected.suffix.lower() in _BLAST_DB_COMPONENT_SUFFIXES:
            selected = selected.with_suffix("")
        self.vars[key].set(str(selected))

    def _start_run(self) -> None:
        try:
            raw_values = {key: variable.get() for key, variable in self.vars.items()}
        except (tk.TclError, TypeError, ValueError) as exc:
            messagebox.showerror("入力エラー", f"入力値を読み取れません。数値欄を確認してください。\n\n{exc}")
            return

        missing = self._missing_required_fields(raw_values)
        if missing:
            messagebox.showerror(
                "AIRR IgBLAST",
                "次の項目を入力してください:\n\n" + "\n".join(f"- {name}" for name in missing),
            )
            return
        try:
            values = _normalize_run_values(raw_values)
        except (TypeError, ValueError) as exc:
            messagebox.showerror("入力エラー", str(exc))
            return

        query_fasta = str(values["query_fasta"]).strip() or None
        try:
            planned_outputs = planned_cpm_output_paths(
                str(values["out"]),
                query_fasta,
            )
        except Exception as exc:
            messagebox.showerror("AIRR IgBLAST", str(exc))
            return
        preflight_errors = _preflight_run_errors(values, planned_outputs)
        if preflight_errors:
            messagebox.showerror(
                "解析開始前の確認",
                "解析を開始できません。次の項目を確認してください:\n\n"
                + "\n".join(f"- {error}" for error in preflight_errors),
            )
            return
        existing = [path for path in planned_outputs if path.exists()]
        overwrite = False
        if existing:
            preview = "\n".join(f"- {path.name}" for path in existing[:8])
            if len(existing) > 8:
                preview += f"\n- ... ほか {len(existing) - 8} ファイル"
            overwrite = messagebox.askyesno(
                "既存結果の置換確認",
                "次の既存結果があります。解析全体が成功した場合だけ置換します。\n\n"
                f"今回のR1: {values['r1']}\n"
                f"今回のR2: {values['r2']}\n\n"
                f"置換対象:\n{preview}\n\n"
                "異なるフォルダの同名sampleを誤って置換していないか確認してください。続行しますか？",
            )
            if not overwrite:
                return
        self.running = True
        self.run_button.configure(state="disabled")
        self._log("Starting IgBLAST run...")
        try:
            thread = threading.Thread(target=self._run_pipeline, args=(values, overwrite), daemon=True)
            thread.start()
        except Exception as exc:
            self.running = False
            self.run_button.configure(state="normal")
            self._log("ERROR: worker thread could not start: " + str(exc))
            messagebox.showerror("AIRR IgBLAST", f"解析スレッドを開始できませんでした。\n\n{exc}")

    def _missing_required_fields(self, values: dict[str, object]) -> list[str]:
        required = [
            ("R1 FASTQ", "r1"),
            ("R2 FASTQ", "r2"),
            ("Output TSV", "out"),
            ("igblastn", "igblastn"),
            ("V DB prefix", "germline_db_v"),
            ("D DB prefix", "germline_db_d"),
            ("J DB prefix", "germline_db_j"),
            ("Aux file", "auxiliary_data"),
        ]
        return [label for label, key in required if not str(values.get(key, "")).strip()]

    def _run_pipeline(self, values: dict[str, object], overwrite: bool) -> None:
        try:
            extra_query_fasta = str(values["query_fasta"]).strip()
            multi_result = run_cpm_umi_igblast_outputs(
                r1_path=str(values["r1"]),
                r2_path=str(values["r2"]),
                output_tsv=str(values["out"]),
                query_fasta=extra_query_fasta or None,
                igblast_config=IgBlastConfig(
                    igblastn=str(values["igblastn"]),
                    germline_db_v=str(values["germline_db_v"]),
                    germline_db_d=str(values["germline_db_d"]) or None,
                    germline_db_j=str(values["germline_db_j"]),
                    auxiliary_data=str(values["auxiliary_data"]) or None,
                    organism=str(values["organism"]),
                    domain_system=str(values["domain_system"]),
                    ig_seqtype=str(values["ig_seqtype"]),
                    num_threads=int(values["num_threads"]),
                ),
                read_selection=str(values["read_selection"]),
                r1_transform=ReadTransform(
                    str(values["r1_orientation"]),
                    int(values["trim_left_r1"]),
                    int(values["trim_right_r1"]),
                ),
                r2_transform=ReadTransform(
                    str(values["r2_orientation"]),
                    int(values["trim_left_r2"]),
                    int(values["trim_right_r2"]),
                ),
                min_length=int(values["min_length"]),
                max_n_rate=float(values["max_n_rate"]),
                query_name_template=str(values["query_name_template"]),
                strict_ids=bool(values["strict_ids"]),
                umi_anchor_max_mismatches=int(values["umi_anchor_max_mismatches"]),
                igblast_batch_size=int(values["igblast_batch_size"]),
                progress_callback=lambda message: self.messages.put(("log", message)),
                work_dir=default_work_dir(),
                overwrite=overwrite,
            )
        except Exception as exc:
            self.messages.put(("error", str(exc)))
            return

        message_parts = []
        for named_result in multi_result.runs:
            result = named_result.result
            stats = result.stats
            pair_stats = result.pair_summary_stats
            message = (
                f"{named_result.label}: {Path(result.output_tsv)}\n"
                f"pairs={stats.total_pairs}, records={stats.records_written}, "
                f"R1={stats.r1_written}, R2={stats.r2_written}, "
                f"skipped_short={stats.skipped_too_short}, skipped_N={stats.skipped_n_rate}, "
                f"UMI_extracted={stats.umi_extracted_pairs}, UMI_missing={stats.umi_missing_pairs}"
            )
            if result.r1_tsv and result.r2_tsv and result.integrated_tsv:
                message += (
                    f"\nR1 TSV: {result.r1_tsv}"
                    f"\nR2 TSV: {result.r2_tsv}"
                    f"\nIntegrated TSV: {result.integrated_tsv}"
                )
            if result.counts_tsv:
                message += f"\nCounts TSV: {result.counts_tsv}"
            if result.counts_xlsx:
                message += f"\nCounts Excel: {result.counts_xlsx}"
            if getattr(result, "umi_counts_tsv", None):
                message += f"\nExact UMI Counts TSV: {result.umi_counts_tsv}"
            if getattr(result, "umi_counts_xlsx", None):
                message += f"\nExact UMI Counts Excel: {result.umi_counts_xlsx}"
            if result.final_productive_counts_tsv:
                message += f"\nFinal productive Counts TSV: {result.final_productive_counts_tsv}"
            if result.final_productive_counts_xlsx:
                message += f"\nFinal productive Counts Excel: {result.final_productive_counts_xlsx}"
            if getattr(result, "final_productive_umi_counts_tsv", None):
                message += (
                    "\nFinal productive exact UMI Counts TSV: "
                    f"{result.final_productive_umi_counts_tsv}"
                )
            if getattr(result, "final_productive_umi_counts_xlsx", None):
                message += (
                    "\nFinal productive exact UMI Counts Excel: "
                    f"{result.final_productive_umi_counts_xlsx}"
                )
            if pair_stats:
                message += (
                    f"\nintegrated_pairs={pair_stats.total_pairs}, "
                    f"junction_aa_conflicts={pair_stats.junction_aa_conflicts}, "
                    f"included_in_counts={pair_stats.included_in_counts}, "
                    f"unique_final_clonotypes={pair_stats.unique_final_clonotypes}, "
                    f"final_productive_included={pair_stats.final_productive_included_in_counts}, "
                    f"unique_final_productive_clonotypes={pair_stats.unique_final_productive_clonotypes}"
                )
            message_parts.append(message)
        message = (
            "Done\n\n"
            f"Completion manifest: {multi_result.manifest_path}\n\n"
            + "\n\n".join(message_parts)
        )
        self.messages.put(("done", message))

    def _poll_messages(self) -> None:
        try:
            while True:
                kind, message = self.messages.get_nowait()
                if kind == "error":
                    self.running = False
                    self._log("ERROR: " + message)
                    messagebox.showerror("AIRR IgBLAST", message)
                    self.run_button.configure(state="normal")
                elif kind == "log":
                    self._log(message)
                elif kind == "done":
                    self.running = False
                    self._log(message)
                    messagebox.showinfo("AIRR IgBLAST", message)
                    self.run_button.configure(state="normal")
        except queue.Empty:
            pass
        self.after(100, self._poll_messages)

    def _on_close(self) -> None:
        if self.running:
            messagebox.showwarning(
                "解析実行中",
                "解析中はウィンドウを閉じられません。完了またはエラー表示までお待ちください。",
            )
            return
        self.master.destroy()

    def _log(self, message: str) -> None:
        self.log.insert("end", message + "\n")
        self.log.see("end")


def main() -> None:
    root = tk.Tk()
    root.title(APP_TITLE)
    root.geometry("900x760")
    App(root)
    root.mainloop()
