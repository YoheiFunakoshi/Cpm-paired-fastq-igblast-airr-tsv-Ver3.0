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
- `final_productive_counts` applies that additional filter and is a secondary
  sensitivity table.

This is a clonotype definition, not whole-BCR nucleotide identity. Synonymous
CDR3 nucleotide differences and SHM differences outside the key can therefore
belong to the same clonotype.

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

## Data-retention invariants

- UMI family counting never removes a read pair from `integrated.tsv`.
- Rows without a countable BCR clonotype remain in `integrated.tsv` with an
  exclusion reason.
- UMI-missing rows remain in pair-level outputs and their valid clonotype.
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
