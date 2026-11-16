# PEDX format specification, version 0.1

PEDX is a tab-separated text format for pedigrees used in teaching and counselling problems. It keeps the first six columns of the PLINK/LINKAGE PED and FAM formats, so a PEDX file can be turned into a `.fam` file without loss of those columns, and adds named columns for the information that NSGC pedigree drawings need (role, death, pregnancy outcomes, twins, adoption, tested genotypes) together with metadata lines describing the diseases and the questions asked about the family.

## 1. File structure

A PEDX file is UTF-8 text with the extension `.pedx`. It has three parts, in this order:

1. metadata lines, each starting with `##`;
2. exactly one column header line, starting with a single `#`;
3. data lines, one individual per line, cells separated by a single tab character.

Empty lines are ignored. A reader must keep the file's line structure so that a file read and written back without changes is byte-identical (see section 6).

```
##format=PEDX/0.1
##family=<ID=FAM01,description="Teaching example: cystic fibrosis and colour blindness in one family">
##disease=<ID=D1,name="Cystic fibrosis",mode=unknown,allele_freq=0.02,penetrance=1.0,mutation_rate=0,pattern=solid>
##disease=<ID=D2,name="Red-green colour blindness",mode=XR,allele_freq=0.08,penetrance=1.0,pattern=hatched>
##union=<p1=II-2,p2=II-3,consanguinity=none,status=married>
##question=<ID=Q1,type=mode,disease=D1>
##question=<ID=Q2,type=risk,individual=II-2,disease=D1,event=carrier>
##question=<ID=Q3,type=risk,individual=III-1,disease=D1,event=affected>
#family	individual	father	mother	sex	D1	D2	G:D1	G:D2	role	deceased	age	kind	order
FAM01	I-1	0	0	1	1	1	.	.	.	1	67	birth	1
FAM01	I-2	0	0	2	1	1	.	.	.	0	64	birth	2
FAM01	II-1	I-1	I-2	1	2	1	1/1	.	P	0	30	birth	1
FAM01	II-2	I-1	I-2	2	1	1	.	.	C	0	27	birth	2
FAM01	II-3	0	0	1	1	2	.	1	.	0	29	birth	3
FAM01	III-1	II-3	II-2	0	0	0	.	.	.	0	.	planned	1
```

## 2. Metadata lines

A metadata line has the form `##key=value` or `##key=<k1=v1,k2=v2,...>`. A structured value may contain commas if it is enclosed in double quotes. Unknown keys, and unknown fields inside a known key, are not errors; they are kept and written back.

| Key | Fields | Notes |
|---|---|---|
| `format` | `PEDX/0.1` | Required; must be the first line. |
| `family` | `ID`, `description` | Optional. |
| `disease` | see section 2.1 | One line per disease. The first disease is the primary disease. |
| `union` | `p1`, `p2`, `consanguinity`, `status` | Only needed for partnerships that cannot be derived from the data lines: childless, separated, or consanguineous with the common ancestor not drawn. `consanguinity` is free text; any value other than `none` marks the union as consanguineous. `status` is `married` (default), `separated` or `childless`. |
| `question` | `ID`, `type`, `individual`, `disease`, `event` | `type` is `mode` (mode-of-inheritance question) or `risk` (risk question for `individual`). `event` is one of `carrier`, `mutant`, `affected`, `intermediate`, `genotype:G`. |

### 2.1. Disease fields

| Field | Values | Default | Used by |
|---|---|---|---|
| `ID` | identifier, also the name of the phenotype column | required | all |
| `name` | text | none | reports, legend |
| `mode` | `AD`, `AR`, `ACD` (autosomal codominant), `XD`, `XR`, `Y`, `MT`, `MF` (multifactorial liability-threshold), `unknown` | `unknown` | Bayes risk, validation (E8, E9, W2), legend |
| `allele_freq` | probability q of the disease allele | 0.01 at computation time | inference, risk |
| `penetrance` | probability that a carrier of the disease genotype is affected | 1.0 | inference, risk, W2 |
| `phenocopy` | probability that a non-carrier is affected | 0 | inference, risk, W2 |
| `mutation_rate` | new mutation rate per transmitted allele (normal → disease) | 0 | inference, risk |
| `prevalence` | population prevalence K | none | `MF` only |
| `heritability` | heritability of liability h², 0 ≤ h² < 1 | none | `MF` only |
| `prevalence_male`, `prevalence_female` | sex-specific prevalence | none | `MF` only |
| `pattern` | `solid`, `hatched`, `dotted`, `quadrant` | `solid` | drawing |

## 3. Column header

The header starts with `#`. The first five columns are positional and identical to PED: `family`, `individual`, `father`, `mother`, `sex`. The sixth column must be the phenotype column of the primary disease. All other columns are recognised by name and may appear in any order. Unknown columns are kept and written back.

## 4. Individual columns

| Column | Values | Missing value |
|---|---|---|
| `family` | family ID | required |
| `individual` | individual ID, unique within the family | required |
| `father`, `mother` | individual ID within the same family, or `0` | `0` |
| `sex` | `1` male, `2` female, `0` unknown | `0` |
| `D1`, `D2`, ... | `1` unaffected, `2` affected (full phenotype), `3` intermediate phenotype (heterozygote expression under codominance; drawn as a half-filled symbol), `0` unknown | `0` |
| `G:D1`, ... | tested genotype: `0/0`, `0/1`, `1/1`; for a male under an X-linked mode, and under `Y` or `MT`, the haploid `0` or `1`. `1` is the disease allele. `1/0` is not accepted. | `.` |
| `role` | `P` proband, `C` consultand | `.` |
| `deceased` | `0`, `1` | `0` |
| `age` | non-negative integer | `.` |
| `kind` | `birth`, `planned`, `pregnancy`, `miscarriage`, `stillbirth` | `birth` |
| `twin` | `MZ:a`, `DZ:a`, `?:a`; the letter names the twin group | `.` |
| `adopted` | `in` (adopted into the family), `out` (adopted out) | `.` |
| `order` | integer, left-to-right order among sibs | `.` |
| `gender` | text, drawn as a label below the symbol (Bennett et al. 2022) | `.` |

## 5. Validation rules

A validator collects every finding instead of stopping at the first one. Each finding has a level, a rule code, the line number (1-based, as shown by an editor), the individual concerned if any, and a message. Level ERROR means the file cannot be used reliably for computation; level WARNING means it can be used but should be checked.

ERROR rules:

| Code | Condition |
|---|---|
| E1 | The first line is not `##format=PEDX/...`, or there is no column header. |
| E2 | A required column is missing, the first five columns are out of order, the sixth column is not the primary disease, or a column name is repeated. |
| E3 | Duplicate individual ID within a family. |
| E4 | `father`, `mother`, `##union` `p1`/`p2` refers to an individual not in the file. |
| E5 | A father has `sex=2`, or a mother has `sex=1`, or the same individual is both father and mother. |
| E6 | The pedigree has a loop (someone is their own ancestor). |
| E7 | A `##disease` has no phenotype column, or a `Dn`/`G:Dn` column has no `##disease` line. |
| E8 | A value outside its allowed set, a malformed metadata line, a line with the wrong number of cells, or a genotype of the wrong ploidy for the disease's mode. |
| E9 | A two-allele genotype for a male under an X-linked mode (`XD`, `XR`). |

WARNING rules:

| Code | Condition |
|---|---|
| W1 | MZ twins differ in sex or in tested genotype. |
| W2 | A tested genotype contradicts the observed phenotype under the disease's mode; checked only for Mendelian modes with `penetrance=1` and `phenocopy=0`. |
| W3 | A `##question` refers to an individual or disease that does not exist. |
| W4 | Only one parent is given (the other is `0`). |

## 6. Writing

A conforming writer writes a file that was read and not changed back byte for byte, including column order, unknown columns, unknown metadata, quoting and the presence or absence of a final newline. Only the lines of changed records are regenerated.

## 7. PLINK `.fam` export

The `.fam` export writes the first six columns: family, individual, father, mother, sex, phenotype. The phenotype comes from a chosen disease (default: the primary disease) and is coded `1` unaffected, `2` affected, `0` missing. The intermediate phenotype `3` is exported as `2`. Missing values are written as `0`, not `-9`. IDs containing whitespace cannot be exported.

## 8. Versioning

The version in the `##format` line follows `PEDX/MAJOR.MINOR`; 0.1 is the only version so far. The pedkit 0.1 reader requires the `PEDX/` prefix (E1 otherwise) and does not check the version number. New optional columns and metadata fields are meant to be added in minor versions; because unknown columns and fields are preserved, older readers keep them intact.
