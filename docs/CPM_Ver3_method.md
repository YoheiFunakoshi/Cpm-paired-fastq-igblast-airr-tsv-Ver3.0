# CPM Ver3.0 analysis method

## Purpose

CPM Ver3.0 separates two observation units that Ver2.0 mixed across two
independent IgBLAST runs:

1. retained R1/R2 read-pair support, comparable to the RG paired-FASTQ tool;
2. exact raw UMI-family support counted within each annotated BCR clonotype.

IgBLAST is run once. Every retained read pair is annotated before UMI family
counting. No R1/R2 full-read Hamming-distance collapse is performed.

## Processing order

```text
paired FASTQ
  -> extract and record the CPM R2 UMI (do not collapse reads)
  -> run R1 and R2 through IgBLAST
  -> integrate the R1/R2 AIRR rows
  -> assign a BCR clonotype
  -> count exact raw UMIs independently inside that clonotype
```

UMI extraction occurs before IgBLAST because the UMI is physically present in
R2. UMI grouping and counting occur after annotation.

## BCR clonotype definition

Ver3.0 uses the same non-merged R1/R2 integration and clonotype definition as
the RG paired-FASTQ tool:

```text
unique_v_gene_set + unique_j_gene_set + final_junction_aa
```

- V/J allele suffixes are removed.
- Candidate sets must match exactly; partial overlap is not a match.
- D calls are retained in the pair-level audit table but are not part of the
  clonotype key.
- `final_junction_aa` must be canonical: present, no stop symbol, starts with
  C, ends with W or F, and is 5-40 amino acids long.
- `integrated_counts` does not require `final_productive=true`.
- `final_productive_counts` applies that additional filter and is a derived
  productive-filtered table.

This is a clonotype definition, not whole-BCR nucleotide identity. Synonymous
CDR3 nucleotide differences and SHM differences outside the key can therefore
belong to the same clonotype.

## Productive annotation and paired-read integration

Ver3.0 does not independently infer productivity. It preserves the
`productive` annotation emitted by IgBLAST for each R1 and R2 AIRR row.
IgBLAST predicts a rearrangement as productive from properties that include an
in-frame V(D)J junction, no internal stop codon, and no internal V-region frame
shift. This is a sequence-based prediction of coding potential, not evidence
that a protein was translated, folded, or expressed on a cell surface.

`final_productive` is neither an R1-and-R2 logical AND nor a logical OR. R1
and R2 are partial observations of the same input molecule. Ver3.0 applies the
following rule:

1. prefer the read selected for `final_junction_aa`;
2. use that read's non-empty productive value;
3. if it is missing, fall back to the opposite read;
4. when both reads have the same junction amino-acid sequence but discordant
   productive values, prefer R1;
5. leave `final_productive` empty when neither read has a value.

Consequently, a pair can have `final_productive=true` when the opposite read is
false or missing. Missing does not mean false: it means that IgBLAST did not
return a productive prediction for that partial read.
AIRR TSV normally represents these values as `T`, `F`, and an empty field;
`true` and `false` in this document denote the corresponding logical values.

`final_junction_aa` is the AIRR JUNCTION amino-acid sequence, including the
conserved C and W/F residues; it is not a full BCR protein sequence. The
productive-filtered workbooks do not additionally require
`complete_vdj=true`.

The raw CPM R2 layout consumes 38 of the 301 sequenced bases with the anchor,
UMI blocks, separators, and `TCTT`, leaving at most about 263 bases of BCR
insert. R2 can therefore receive a V call without reaching enough J/junction
sequence for a productive prediction. Requiring both reads to be productive is
an optional stringent QC, not the Ver3.0 biological definition.

References: [AIRR Rearrangement Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
and [NCBI IgBLAST productive update](https://blast.ncbi.nlm.nih.gov/doc/blast-news/2021-BLAST-News.html).

## Exact-UMI rule

Within each clonotype, each distinct observed raw 12-mer containing only A, C,
G, and T contributes one UMI family. Repeated read pairs with the same UMI and
the same clonotype contribute one family but remain present in the pair-level
output. Ambiguous or malformed UMIs are retained as UMI-missing read-pair
support rather than being called molecular families.

The same UMI text occurring in another clonotype is processed independently.
There is no sample-wide UMI merge and no shared-UMI decision rule.

Ver3.0 does not merge one-base-neighbor UMIs. It records the raw UMI so that a
future, separately named error-correction mode can be evaluated without
changing the Ver3.0 default semantics.

## UMI unavailable

If a usable UMI cannot be extracted, the read pair is not discarded. It
contributes to `umi_missing_read_pair_count` for its clonotype. Because PCR
duplicates cannot be resolved without a UMI, this is a read-pair support unit,
not a molecular-family unit.

The convenience value is therefore named `inclusive_support_count`:

```text
inclusive_support_count
  = umi_family_count + umi_missing_read_pair_count
```

It must not be described as a strict molecule count. The two components and
`read_pair_count` are always reported separately.

### Two different UMI-missing percentages

The extraction and workbook metrics come from different processing stages and
use different denominators:

```text
input-pair UMI-assignment failure rate (prepare stage)
  = umi_missing_pairs / total_pairs

workbook pair-level UMI-missing rate
  = sum(umi_missing_read_pair_count) / sum(read_pair_count)

UMI-missing share of inclusive support
  = sum(umi_missing_read_pair_count) / sum(inclusive_support_count)
```

In one validation run, 96,030 of 644,824 input pairs (14.892%) had no usable
UMI assignment. In the productive subset, 64,178 of 413,106 countable pairs
(15.535%) were UMI-missing. The 348,928 UMI-known pairs were reduced to 74,968
clonotype-local exact UMI families, whereas all 64,178 UMI-missing pairs
remained one-for-one:

```text
inclusive_support_count = 74,968 + 64,178 = 139,146
UMI-missing share of inclusive support = 64,178 / 139,146 = 46.123%
```

Thus, 46.123% is not the UMI extraction-failure rate. It is a denominator
effect caused by collapsing only the UMI-known side (approximately 4.65 known
pairs per family) while retaining missing pairs without molecular correction.

In a separate audit of the raw FASTQ, all R2 reads in this validation run were
301 nt. No anchor-passing read had an ambiguous 12-mer; every missing
assignment resulted from the expected anchor
exceeding the allowed two mismatches, and 89,481 of 96,030 missing pairs
(93.18%) differed at six or more anchor positions. This is reported as failure
to confirm the expected anchor, not as proof that the 12 UMI bases themselves
were unreadable. Ordinary one- or two-base sequencing errors alone are
unlikely to explain most missing assignments; alternative library structures,
positional shifts, and off-target products are possible, but this FASTQ alone
does not establish the cause.

There is no assay-independent normal UMI-extraction rate. Approximately 85%
usable UMI pairs in this run is not treated as an automatic run failure, but
the pair-level missing rate must be compared across samples produced with the
same protocol. Large between-sample differences can bias
`inclusive_support_count` and require review of library structure, anchor
settings, read orientation, and sequencing quality. MIGEC similarly recommends
investigating the cause of low barcode extraction. The BCR-oriented abstar
workflow likewise leaves the UMI
empty but continues normal annotation when a conserved pattern cannot be
matched. A separate 10x Genomics single-cell assay lists >75% valid UMIs as
ideal, but that assay-specific value is context only and is not used as a CPM
BCR acceptance threshold.

References: [abstar UMI support](https://abstar.readthedocs.io/en/stable/umis.html),
[MIGEC documentation](https://migec.readthedocs.io/_/downloads/en/latest/pdf/),
and [10x Genomics assay-specific QC example](https://assets.ctfassets.net/an68im79xiti/1R7z9gj36IkuqRpdo2W8Un/fe76580732b3de5c449340278975a658/CG000475_TechNote_ChromiumNextGEM_SC3-_CMOWebSummary_Rev_B.pdf).

## Strict exact-UMI-family view

The inclusive UMI workbooks remain unchanged and retain UMI-missing support.
The same run also writes two convenience views for analyses that require a
denominator made only from valid exact UMI families:

- `exact_umi_family_counts.xlsx` without a productive restriction, including
  true, false, and missing productive annotations where otherwise countable;
- `final_productive_exact_umi_family_counts.xlsx` for the independently
  recalculated productive subset.

Only clonotypes with `umi_family_count > 0` are included. Each strict workbook
contains five columns: normalized V and J candidate sets, canonical final
junction amino acid, `umi_family_count`, and `umi_family_percent`.

```text
umi_family_percent
  = umi_family_count / sum(umi_family_count in that strict workbook)
```

UMI-missing pairs contribute neither to the primary count nor to this
denominator. A clonotype supported only by UMI-missing pairs is therefore absent
from the strict view. This is an additional summary view, not deletion from
`integrated.tsv`, and it does not change either inclusive UMI workbook.

## Examples

For one clonotype:

```text
UMI A, UMI A, UMI A, UMI B, missing, missing
```

the results are:

```text
read_pair_count = 6
umi_family_count = 2
umi_missing_read_pair_count = 2
inclusive_support_count = 4
```

For another clonotype:

```text
UMI C, UMI C, UMI D, UMI E, UMI A, UMI A, missing
```

the results are:

```text
read_pair_count = 7
umi_family_count = 4
umi_missing_read_pair_count = 1
inclusive_support_count = 5
```

UMI A is counted once in each clonotype. Neither occurrence is discarded or
globally reconciled.

## Six Excel workbooks

The six workbooks are not independent IgBLAST runs. They show read-pair,
inclusive UMI, and strict exact-UMI-family views, each with and without the
productive restriction.

| Workbook | Main support counts reported | Included pairs | Primary use |
|---|---|---|---|
| `integrated_counts.xlsx` | read pair | V, J, and canonical junction AA present | RG-compatible read-pair comparison and audit |
| `final_productive_counts.xlsx` | read pair | basic counts plus `final_productive=true` | productive-filtered RG comparison |
| `umi_counts.xlsx` | exact raw UMI family, UMI-missing pair, and read pair | basic counts | CPM UMI-family support and PCR-duplication context |
| `final_productive_umi_counts.xlsx` | exact raw UMI family, UMI-missing pair, and read pair | `final_productive=true` | productive-filtered CPM repertoire evaluation |
| `exact_umi_family_counts.xlsx` | exact raw UMI family | `umi_family_count > 0` | strict valid-UMI-family abundance and percentage |
| `final_productive_exact_umi_family_counts.xlsx` | exact raw UMI family | productive subset with `umi_family_count > 0` | productive-filtered strict UMI-family evaluation |

Each row is one clonotype key:

```text
unique_v_gene_set + unique_j_gene_set + final_junction_aa
```

The row unit is a clonotype, while `read_pair_count`, `umi_family_count`, and
`umi_missing_read_pair_count` are distinct support units.

The two read-pair workbooks contain the same ten columns:

- normalized V and J candidate sets and canonical final junction amino acid;
- `read_pair_count`;
- R1/R2 integration-status counts (`match`, `conflict`, `r1_only`, and
  `r2_only`);
- `productive_true_count` and `canonical_junction_aa_count`.

`read_pair_count` can contain PCR duplicates and is not a molecule count.

The two inclusive UMI workbooks contain 15 columns. In addition to the clonotype and
integration-status fields, they report:

- `umi_family_count`: distinct valid exact raw 12-mer UMIs in the clonotype;
- `read_pair_count`: all retained countable pairs in the clonotype;
- `umi_known_read_pair_count`: pairs with a valid UMI, including PCR copies;
- `umi_missing_read_pair_count`: retained pairs without a usable UMI;
- `inclusive_support_count`: UMI families plus UMI-missing pairs;
- `inclusive_support_percent`: the share of the workbook's total inclusive
  support.

`inclusive_support_count` is deliberately a hybrid support measure, not a
strict molecule count. `final_productive_umi_counts.xlsx` independently
recomputes its UMI sets, missing counts, inclusive support, and percentages
from productive pairs; it is not merely a row filter applied to
`umi_counts.xlsx`.

The two strict exact-UMI-family workbooks contain only five columns:

- normalized V and J candidate sets and canonical final junction amino acid;
- `umi_family_count`;
- `umi_family_percent`, whose denominator is the sum of `umi_family_count` in
  that strict workbook.

They contain only `umi_family_count > 0` rows and are sorted primarily by
`umi_family_count` descending. Their TSV companions use the same five-column
schema. `final_productive_exact_umi_family_counts` is derived from UMI families
recomputed in the productive subset before the strict row filter is applied.

This is not a productive consensus call for an entire UMI family. If the same
clonotype and UMI contain both productive and non-productive pairs, that UMI is
counted once in the productive subset whenever at least one productive pair is
present.

For a study focused on potentially protein-coding repertoires,
`final_productive_umi_counts.xlsx` can be the primary CPM analysis candidate,
while the unfiltered workbooks remain the audit and sensitivity reference. To
measure the effect of the productive filter, compare like with like:

```text
clonotype retention = productive workbook rows / unfiltered workbook rows
read-pair retention = sum(productive read_pair_count) / sum(unfiltered read_pair_count)
UMI-family retention = sum(productive umi_family_count) / sum(unfiltered umi_family_count)
inclusive retention = sum(productive inclusive_support_count) / sum(unfiltered inclusive_support_count)
```

These denominators are different observation units and must not be
interchanged. Removing low-support singleton clonotypes can reduce clonotype
richness more than it reduces read-pair or UMI-family support.

For QASAS or another matching analysis, an inclusive-table matched-clonotype
count and a strict-table UMI-family percentage may be reported side by side
only as explicitly named, different metrics. The former answers how many
clonotypes were detected while retaining missing-only evidence; the latter
answers what share of valid exact UMI families matched. They must not be
presented as one common observation unit. A strict abundance calculation uses
`umi_family_count` or `umi_family_percent`, never `read_pair_count` or
`inclusive_support_count`.

## Data-retention invariants

- UMI family counting never removes a read pair from `integrated.tsv`.
- Rows without a countable BCR clonotype remain in `integrated.tsv` with an
  exclusion reason.
- UMI-missing rows remain in pair-level outputs and their valid clonotype.
- Strict workbooks omit missing-only clonotypes only from that summary view;
  they do not remove or rewrite pair-level evidence.
- The RG-compatible read-pair tables and UMI tables are produced from the same
  integrated rows, so their provenance cannot diverge through separate
  IgBLAST runs.
- A run is complete only when all outputs and the completion manifest have
  been published successfully.

## Interpretation boundary

This method is annotation-first, clonotype-aware UMI counting. It is not a UMI
consensus caller and does not construct an error-corrected molecule sequence.
An annotation error that changes V, J, or junction amino acid can split one
true molecule across clonotypes. The pair-level output is retained so that
such cases can be audited, and a future consensus mode must be introduced as a
new, explicitly versioned method rather than silently changing these counts.
