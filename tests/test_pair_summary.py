from __future__ import annotations

import csv
from contextlib import contextmanager
from pathlib import Path
import shutil
import unittest
import uuid
import zipfile

from airr_igblast_paired.pair_summary import (
    COUNTS_FIELDNAMES,
    UMI_COUNTS_FIELDNAMES,
    _format_percent,
    default_derived_tsv_paths,
    extract_umi,
    pair_id_and_read_label,
    split_and_integrate_airr_tsv,
)


AIRR_HEADER = (
    "sequence_id\tv_call\td_call\tj_call\tjunction\tjunction_aa\tproductive"
)
RG_COUNTS_FIELDNAMES = [
    "unique_v_gene_set",
    "unique_j_gene_set",
    "final_junction_aa",
    "read_pair_count",
    "match_count",
    "conflict_count",
    "r1_only_count",
    "r2_only_count",
    "productive_true_count",
    "canonical_junction_aa_count",
]


def write_airr(path: Path, rows: list[str]) -> None:
    path.write_text("\n".join([AIRR_HEADER, *rows]) + "\n", encoding="utf-8")


def read_tsv(path: Path) -> tuple[list[str] | None, list[dict[str, str]]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


@contextmanager
def test_directory():
    root = Path(f"test_pair_summary_tmp_{uuid.uuid4().hex[:8]}")
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir()
    try:
        yield root
    finally:
        shutil.rmtree(root, ignore_errors=True)


class PairSummaryTests(unittest.TestCase):
    def test_airr_header_must_include_required_columns(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "broken.airr.tsv"
            input_tsv.write_text("sequence_id\tv_call\nread1|R1\tIGHV1\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "missing required column"):
                split_and_integrate_airr_tsv(input_tsv)

    def test_airr_rows_must_follow_query_name_contract(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "broken.airr.tsv"
            write_airr(input_tsv, ["read1\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCARYW\tT"])

            with self.assertRaisesRegex(ValueError, "R1/R2 query-name contract"):
                split_and_integrate_airr_tsv(input_tsv)

    def test_airr_duplicate_read_side_is_rejected(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "broken.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "read1|R1\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCARYW\tT",
                    "read1|R1\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCARYW\tT",
                ],
            )

            with self.assertRaisesRegex(ValueError, "duplicate R1 row"):
                split_and_integrate_airr_tsv(input_tsv)

    def test_default_derived_paths_include_rg_and_umi_outputs(self) -> None:
        paths = default_derived_tsv_paths(Path("results") / "sample.airr.tsv")

        self.assertEqual(paths.r1_tsv, Path("results") / "sample.R1.airr.tsv")
        self.assertEqual(paths.r2_tsv, Path("results") / "sample.R2.airr.tsv")
        self.assertEqual(paths.integrated_tsv, Path("results") / "sample.integrated.tsv")
        self.assertEqual(paths.counts_tsv, Path("results") / "sample.integrated_counts.tsv")
        self.assertEqual(paths.counts_xlsx, Path("results") / "sample.integrated_counts.xlsx")
        self.assertEqual(
            paths.final_productive_counts_tsv,
            Path("results") / "sample.final_productive_counts.tsv",
        )
        self.assertEqual(
            paths.final_productive_counts_xlsx,
            Path("results") / "sample.final_productive_counts.xlsx",
        )
        self.assertEqual(paths.umi_counts_tsv, Path("results") / "sample.umi_counts.tsv")
        self.assertEqual(paths.umi_counts_xlsx, Path("results") / "sample.umi_counts.xlsx")
        self.assertEqual(
            paths.final_productive_umi_counts_tsv,
            Path("results") / "sample.final_productive_umi_counts.tsv",
        )
        self.assertEqual(
            paths.final_productive_umi_counts_xlsx,
            Path("results") / "sample.final_productive_umi_counts.xlsx",
        )

    def test_pair_id_and_umi_parsing_follow_cpm_query_name_contract(self) -> None:
        self.assertEqual(pair_id_and_read_label("read-1|R1"), ("read-1", "R1"))
        self.assertEqual(
            pair_id_and_read_label("read-1|R2|UMI=ACGTACGTACGT"),
            ("read-1|UMI=ACGTACGTACGT", "R2"),
        )
        self.assertEqual(pair_id_and_read_label("read-1"), ("read-1", None))
        self.assertEqual(extract_umi("read-1|R2|UMI=ACGTACGTACGT"), "ACGTACGTACGT")
        self.assertEqual(
            extract_umi("UMI=READID|R2|UMI=ACGTACGTACGT"),
            "ACGTACGTACGT",
        )
        self.assertEqual(extract_umi("read-1|R2|UMI=NA"), "")

    def test_rg_count_outputs_keep_rg_exact_ten_column_schema(self) -> None:
        self.assertEqual(COUNTS_FIELDNAMES, RG_COUNTS_FIELDNAMES)
        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "read1|R1|UMI=AAAAAAAAAAAA\tIGHV1*01\tIGHD1\tIGHJ4*01\tAAA\tCARYW\tT",
                    "read1|R2|UMI=AAAAAAAAAAAA\tIGHV1*02\tIGHD1\tIGHJ4*02\tAAA\tCARYW\tT",
                    "read2|R1|UMI=AAAAAAAAAAAA\tIGHV1*03\tIGHD1\tIGHJ4*03\tAAA\tCARYW\tT",
                    "read2|R2|UMI=AAAAAAAAAAAA\tIGHV1*04\tIGHD1\tIGHJ4*04\tAAA\tCARYW\tT",
                ],
            )

            paths, _stats = split_and_integrate_airr_tsv(input_tsv)
            fields, rows = read_tsv(paths.counts_tsv)
            productive_fields, productive_rows = read_tsv(paths.final_productive_counts_tsv)

            self.assertEqual(fields, RG_COUNTS_FIELDNAMES)
            self.assertEqual(productive_fields, RG_COUNTS_FIELDNAMES)
            self.assertEqual(rows[0]["read_pair_count"], "2")
            self.assertEqual(productive_rows[0]["read_pair_count"], "2")
            self.assertNotIn("umi_family_count", rows[0])

    def test_umi_counting_is_independent_inside_each_bcr_pattern(self) -> None:
        """The same raw UMI is one family in each BCR key, never globally."""

        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            rows: list[str] = []

            # User example, pattern 1: A x3, B x1, missing x2 -> 2 + 2 = 4.
            for index in range(3):
                rows.append(
                    f"p1_a{index}|R2|UMI=AAAAAAAAAAAA\tIGHV1*01\tIGHD1\tIGHJ4*01\tAAA\tCPATTERNONEW\tT"
                )
            rows.extend(
                [
                    "p1_b|R2|UMI=ACGTACGTACGT\tIGHV1*02\tIGHD1\tIGHJ4*02\tAAA\tCPATTERNONEW\tT",
                    "p1_missing1|R2|UMI=NA\tIGHV1*01\tIGHD1\tIGHJ4*01\tAAA\tCPATTERNONEW\tT",
                    "p1_missing2|R2|UMI=NA\tIGHV1*01\tIGHD1\tIGHJ4*01\tAAA\tCPATTERNONEW\tT",
                ]
            )

            # Pattern 2: C x2, D, E, and the same A x2, plus missing -> 4 + 1 = 5.
            for index in range(2):
                rows.append(
                    f"p2_c{index}|R2|UMI=CCCCCCCCCCCC\tIGHV2*01\tIGHD2\tIGHJ6*01\tCCC\tCPATTERNTWOF\tT"
                )
                rows.append(
                    f"p2_a{index}|R2|UMI=AAAAAAAAAAAA\tIGHV2*01\tIGHD2\tIGHJ6*01\tCCC\tCPATTERNTWOF\tT"
                )
            rows.extend(
                [
                    "p2_d|R2|UMI=GGGGGGGGGGGG\tIGHV2*01\tIGHD2\tIGHJ6*01\tCCC\tCPATTERNTWOF\tT",
                    "p2_e|R2|UMI=TTTTTTTTTTTT\tIGHV2*01\tIGHD2\tIGHJ6*01\tCCC\tCPATTERNTWOF\tT",
                    "p2_missing|R2|UMI=NA\tIGHV2*01\tIGHD2\tIGHJ6*01\tCCC\tCPATTERNTWOF\tT",
                ]
            )
            # An unassignable BCR row remains in integrated.tsv, not in count tables.
            rows.append("unassigned|R2|UMI=AAAAAAAAAAAA\t\t\t\t\t\tF")
            write_airr(input_tsv, rows)

            paths, stats = split_and_integrate_airr_tsv(input_tsv)
            _fields, integrated = read_tsv(paths.integrated_tsv)
            count_fields, read_pair_rows = read_tsv(paths.counts_tsv)
            umi_fields, umi_rows = read_tsv(paths.umi_counts_tsv)

            self.assertEqual(stats.total_pairs, 14)
            self.assertEqual(len(integrated), 14)
            unassigned = next(row for row in integrated if row["pair_id"].startswith("unassigned"))
            self.assertEqual(unassigned["include_in_counts"], "false")
            self.assertEqual(unassigned["umi"], "AAAAAAAAAAAA")
            self.assertEqual(count_fields, RG_COUNTS_FIELDNAMES)
            self.assertEqual(umi_fields, UMI_COUNTS_FIELDNAMES)

            read_lookup = {row["final_junction_aa"]: row for row in read_pair_rows}
            umi_lookup = {row["final_junction_aa"]: row for row in umi_rows}
            self.assertEqual(read_lookup["CPATTERNONEW"]["read_pair_count"], "6")
            self.assertEqual(read_lookup["CPATTERNTWOF"]["read_pair_count"], "7")

            pattern1 = umi_lookup["CPATTERNONEW"]
            self.assertEqual(pattern1["umi_family_count"], "2")
            self.assertEqual(pattern1["umi_known_read_pair_count"], "4")
            self.assertEqual(pattern1["umi_missing_read_pair_count"], "2")
            self.assertEqual(pattern1["inclusive_support_count"], "4")

            pattern2 = umi_lookup["CPATTERNTWOF"]
            self.assertEqual(pattern2["umi_family_count"], "4")
            self.assertEqual(pattern2["umi_known_read_pair_count"], "6")
            self.assertEqual(pattern2["umi_missing_read_pair_count"], "1")
            self.assertEqual(pattern2["inclusive_support_count"], "5")
            # AAAAAAAAAAAA contributes once to each pattern without shared-UMI logic.
            self.assertEqual(sum(int(row["umi_family_count"]) for row in umi_rows), 6)

    def test_exact_raw_umi_values_are_not_hamming_corrected(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "one|R2|UMI=AAAAAAAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCEXACTW\tT",
                    "two|R2|UMI=AAAAAAAAAAAC\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCEXACTW\tT",
                    "three|R2|UMI=AAAAAAAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCEXACTW\tT",
                    "short|R2|UMI=AAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCEXACTW\tT",
                    "ambiguous|R2|UMI=AAAAANAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCEXACTW\tT",
                ],
            )

            paths, _stats = split_and_integrate_airr_tsv(input_tsv)
            _fields, rows = read_tsv(paths.umi_counts_tsv)

            self.assertEqual(rows[0]["umi_family_count"], "2")
            self.assertEqual(rows[0]["umi_known_read_pair_count"], "3")
            self.assertEqual(rows[0]["umi_missing_read_pair_count"], "2")
            self.assertEqual(rows[0]["inclusive_support_count"], "4")

    def test_final_productive_outputs_recalculate_counts_and_percentages(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "p1_a_t|R2|UMI=AAAAAAAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCPRODONEW\tT",
                    "p1_a_f|R2|UMI=AAAAAAAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCPRODONEW\tF",
                    "p1_b_f|R2|UMI=ACGTACGTACGT\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCPRODONEW\tF",
                    "p1_missing_t|R2|UMI=NA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCPRODONEW\tT",
                    "p2_c_t|R2|UMI=CCCCCCCCCCCC\tIGHV2\tIGHD2\tIGHJ6\tCCC\tCPRODTWOF\tT",
                ],
            )

            paths, _stats = split_and_integrate_airr_tsv(input_tsv)
            _fields, all_read_rows = read_tsv(paths.counts_tsv)
            _fields, productive_read_rows = read_tsv(paths.final_productive_counts_tsv)
            _fields, all_umi_rows = read_tsv(paths.umi_counts_tsv)
            _fields, productive_umi_rows = read_tsv(paths.final_productive_umi_counts_tsv)

            all_read = {row["final_junction_aa"]: row for row in all_read_rows}
            productive_read = {row["final_junction_aa"]: row for row in productive_read_rows}
            all_umi = {row["final_junction_aa"]: row for row in all_umi_rows}
            productive_umi = {row["final_junction_aa"]: row for row in productive_umi_rows}

            self.assertEqual(all_read["CPRODONEW"]["read_pair_count"], "4")
            self.assertEqual(productive_read["CPRODONEW"]["read_pair_count"], "2")
            self.assertEqual(all_umi["CPRODONEW"]["umi_family_count"], "2")
            self.assertEqual(all_umi["CPRODONEW"]["inclusive_support_count"], "3")
            self.assertEqual(all_umi["CPRODONEW"]["inclusive_support_percent"], "75.000000%")
            self.assertEqual(productive_umi["CPRODONEW"]["umi_family_count"], "1")
            self.assertEqual(productive_umi["CPRODONEW"]["umi_known_read_pair_count"], "1")
            self.assertEqual(productive_umi["CPRODONEW"]["umi_missing_read_pair_count"], "1")
            self.assertEqual(productive_umi["CPRODONEW"]["inclusive_support_count"], "2")
            self.assertEqual(
                productive_umi["CPRODONEW"]["inclusive_support_percent"],
                "66.666667%",
            )

    def test_gene_sets_canonical_filter_and_r1_r2_rule_match_rg(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "set_ab|R1|UMI=AAAAAAAAAAAA\tIGHV4-61*01,IGHV4-59*01\tIGHD1\tIGHJ4*01\tAAA\tCVQGFDYW\tT",
                    "set_ab|R2|UMI=AAAAAAAAAAAA\tIGHV4-59*02,IGHV4-61*03\tIGHD2\tIGHJ4*02\tCCC\tCVQGFDYW\tF",
                    "bad_start|R2|UMI=ACGTACGTACGT\tIGHV4-61*01\tIGHD1\tIGHJ4*01\tAAA\tAVQGFDYW\tT",
                    "missing_j|R2|UMI=CCCCCCCCCCCC\tIGHV4-61*01\tIGHD1\t\tAAA\tCVQGFDYW\tT",
                    "canonical_r2|R1|UMI=GGGGGGGGGGGG\tIGHV1\tIGHD1\tIGHJ4\tAAA\tARYFDYW\tT",
                    "canonical_r2|R2|UMI=GGGGGGGGGGGG\tIGHV2\tIGHD1\tIGHJ4\tAAA\tCARYFDYW\tF",
                ],
            )

            paths, stats = split_and_integrate_airr_tsv(input_tsv)
            _fields, integrated_rows = read_tsv(paths.integrated_tsv)
            integrated = {row["pair_id"].split("|UMI=", 1)[0]: row for row in integrated_rows}
            _fields, counts = read_tsv(paths.counts_tsv)

            self.assertEqual(stats.total_pairs, 4)
            self.assertEqual(stats.included_in_counts, 2)
            self.assertEqual(integrated["set_ab"]["unique_v_gene_set"], "IGHV4-59,IGHV4-61")
            self.assertEqual(integrated["set_ab"]["unique_j_gene_set"], "IGHJ4")
            self.assertEqual(integrated["set_ab"]["final_productive"], "T")
            self.assertEqual(integrated["bad_start"]["include_in_counts"], "false")
            self.assertIn("junction_aa_not_c_start", integrated["bad_start"]["exclude_reason"])
            self.assertEqual(integrated["missing_j"]["include_in_counts"], "false")
            self.assertIn("missing_j_call", integrated["missing_j"]["exclude_reason"])
            self.assertEqual(integrated["canonical_r2"]["final_junction_aa"], "CARYFDYW")
            self.assertEqual(integrated["canonical_r2"]["preferred_read"], "R2")
            self.assertEqual(integrated["canonical_r2"]["final_productive"], "F")
            self.assertEqual(len(counts), 2)

    def test_umi_xlsx_percent_is_numeric_and_all_outputs_exist(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "sample.airr.tsv"
            write_airr(
                input_tsv,
                [
                    "known|R2|UMI=AAAAAAAAAAAA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCXLSXW\tT",
                    "missing|R2|UMI=NA\tIGHV1\tIGHD1\tIGHJ4\tAAA\tCXLSXW\tT",
                ],
            )

            paths, _stats = split_and_integrate_airr_tsv(input_tsv)
            for path in (
                paths.counts_xlsx,
                paths.final_productive_counts_xlsx,
                paths.umi_counts_xlsx,
                paths.final_productive_umi_counts_xlsx,
            ):
                self.assertTrue(path.exists())

            with zipfile.ZipFile(paths.umi_counts_xlsx) as archive:
                sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
                styles_xml = archive.read("xl/styles.xml").decode("utf-8")
                app_xml = archive.read("docProps/app.xml").decode("utf-8")
                core_xml = archive.read("docProps/core.xml").decode("utf-8")
                self.assertIn("final_junction_aa (canonical)", sheet_xml)
                self.assertIn('<c r="I2" s="1"><v>1.000000</v></c>', sheet_xml)
                self.assertNotIn('<c r="I2" t="inlineStr"', sheet_xml)
                self.assertIn('formatCode="0.000000%"', styles_xml)
                self.assertIn("Ver3.0", app_xml)
                self.assertIn("Ver3.0", core_xml)
                self.assertNotIn("Ver2.0", app_xml + core_xml)

    def test_empty_input_writes_headers_for_every_count_output(self) -> None:
        with test_directory() as directory:
            input_tsv = Path(directory) / "empty.airr.tsv"
            input_tsv.write_text("", encoding="utf-8")

            paths, stats = split_and_integrate_airr_tsv(input_tsv)

            self.assertEqual(stats.total_pairs, 0)
            for path, expected_fields in (
                (paths.counts_tsv, COUNTS_FIELDNAMES),
                (paths.final_productive_counts_tsv, COUNTS_FIELDNAMES),
                (paths.umi_counts_tsv, UMI_COUNTS_FIELDNAMES),
                (paths.final_productive_umi_counts_tsv, UMI_COUNTS_FIELDNAMES),
            ):
                fields, rows = read_tsv(path)
                self.assertEqual(fields, expected_fields)
                self.assertEqual(rows, [])

    def test_percent_precision_preserves_large_singleton_table(self) -> None:
        singleton_count = 250_000
        singleton_percent = _format_percent(1, singleton_count)

        self.assertEqual(singleton_percent, "0.000400%")
        self.assertNotEqual(singleton_percent, "0.000000%")
        total_percent = float(singleton_percent.removesuffix("%")) * singleton_count
        self.assertAlmostEqual(total_percent, 100.0, places=6)


if __name__ == "__main__":
    unittest.main()
