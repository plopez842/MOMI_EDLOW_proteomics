# Local data inputs

Participant-level data are intentionally excluded from Git. Download the study data from ImmPort accession **SDY2913**, then create the two analysis-ready R objects below.

## `comb_v0_vax.rds`

A data frame with one row per plasma sample. Required metadata columns:

| Column | Meaning |
|---|---|
| `Sample_Label` | Unique sample identifier |
| `Timepoint_v1` | `V0`, `V1`, or `V2` |
| `GA_collection_days` | Gestational age at sample collection, in days |
| `Trim_Collec` | `Trim_1st`, `Trim_2nd`, or `Trim_3rd` |
| `Vax1_Trim` | Trimester of first vaccine dose |
| `diff_collec_vax1` | Days from Dose 1 to collection |
| `diff_collec_vax2` | Days from Dose 2 to collection |

Protein columns must be named with the SOMAmer identifiers used in the study, beginning with `seq.`. Values are log2 RFU.

## `aptamer_annotations.rds`

A data frame with these required columns:

| Column | Meaning |
|---|---|
| `AptName` | SOMAmer identifier matching a protein column |
| `EntrezGeneSymbol` | Display gene symbol |
| `Single_UniPro` | Single UniProt identifier used for pathway analysis |

Optional columns such as `F_stat` can be included. When several aptamers map to the same UniProt identifier, the workflow keeps the aptamer with the largest `F_stat`; if that column is unavailable, the aptamer with the largest variance in the study data is used.

## Placement

```text
data/processed/comb_v0_vax.rds
data/processed/aptamer_annotations.rds
```

The files can live elsewhere if `MOMI_DATA_DIR` points to their directory.
