# pedkit user guide

pedkit reads a pedigree written in the PEDX format, checks it for errors, draws it with NSGC symbols, evaluates the mode of inheritance of a tracked disease by exact computation or Monte Carlo simulation, computes the Bayes risk of a chosen individual, and for multifactorial diseases computes the risk under the liability-threshold model and plots it. This guide explains what each command does, what it prints and how to read the output. Every interpretation section states what a result means and also what it does not mean, because misuse of tools like these usually comes from interpretation, not from the computation.

The full technical specification of the format is in `docs/format.md`. This guide covers what a user needs to know and does not repeat the specification. How the numerical results were checked is described in `docs/validation.md`.

## 1. Installation

pedkit needs Python 3.10 or newer. Reading, writing, validation and SVG drawing use only the standard library. The other features need optional packages, installed as extras:

```sh
pip install pedkit                  # read, write, validate, draw SVG
pip install "pedkit[inference]"     # + mode inference and Bayes risk (numpy)
pip install "pedkit[mf]"            # + multifactorial risk and plots (numpy, matplotlib)
pip install "pedkit[qt]"            # + viewer window, PDF/PNG export (PyQt5)
pip install "pedkit[all]"           # everything
```

From a source checkout, install in editable mode with the test dependencies:

```sh
git clone https://github.com/mkamil-turan/pedkit.git
cd pedkit
python -m venv .venv
.venv/bin/pip install -e ".[all,test]"
```

PDF/PNG export and the viewer render the SVG with Qt's SVG module. On some Linux distributions this module is a separate system package that pip does not install (on Arch Linux `qt5-svg`, on Debian/Ubuntu `libqt5svg5`). If the viewer reports `Could not load PyQt5/QtSvg`, install that package.

The fastest check that the installation works is to run the tests. `python -m pytest -q` should end with a line such as `90 passed` and no `failed`. If a test fails, first check for a missing package.

The commands below are written as `pedkit ...`. The same commands work as `python -m pedkit ...`. The directory layout of the source tree is: `pedkit/` is the package, `examples/` holds the `.pedx` files used throughout this guide, `tests/` holds the automated tests, `docs/` holds the documentation.

## 2. Writing a PEDX file

A PEDX file is UTF-8 plain text with the extension `.pedx`. It has three parts: metadata lines starting with `##`, a single column header line starting with `#`, and data lines that each describe one individual. The columns of a data line are separated by a tab character, not spaces. If your editor converts tabs to spaces the reader cannot split the columns, and validation reports errors such as "Line has 1 cells, the header has 14 columns."

The following file is `examples/cf_family.pedx` and is used for most of this guide (tabs look like spaces here):

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

The first line must always be `##format=PEDX/0.1`. Each `##disease` line defines a disease, and the first disease in the file is the primary disease. When a command is not told which disease to use, it uses the primary disease. The fields of the disease line enter the computations directly, and we come back to them when interpreting inference results. `allele_freq` is the population frequency q of the disease allele; if omitted, 0.01 is assumed. `penetrance` is the probability that a carrier of the disease genotype is affected; if omitted, 1.0. `phenocopy` is the probability that someone without the disease genotype appears affected anyway, and `mutation_rate` is the new mutation rate per generation; both default to 0. The `mode` field can be `AD`, `AR`, `ACD`, `XD`, `XR`, `Y`, `MT`, `MF` or `unknown`. It does not affect mode inference, because inference tries every Mendelian mode. In the Bayes risk only the given mode is used if one is set. `pattern` sets how an affected symbol is filled in the drawing: `solid`, `hatched`, `dotted` or `quadrant`.

The first five columns of the header are fixed and cannot be reordered: `family`, `individual`, `father`, `mother`, `sex`. The sixth column must be the phenotype column of the primary disease, named after its ID (`D1`). The remaining columns are recognised by name and may appear in any order. Unrecognised columns (for example a `note` column of your own) cause no error; they are kept and written back unchanged.

Cell values are written as follows. The `father` and `mother` columns contain the parent's individual ID, or `0` if the parent is not in the file. `sex` is 1 male, 2 female, 0 unknown. In phenotype columns 1 is unaffected, 2 affected, 0 unknown and 3 an intermediate phenotype. The value 3 is meaningful only for codominant diseases and records the heterozygote's own expression; sickle cell trait is an example. Genotype columns carry the disease ID with the prefix `G:` (`G:D1`). Their values are tested genotypes: `0/0`, `0/1`, `1/1`, for a male in an X-linked disease `0` or `1`, and `.` if not tested. Here 1 is the disease allele and 0 the normal allele. The spelling `1/0` is not accepted; a heterozygote is always `0/1`. The `role` column is `P` for the proband, `C` for the consultand, `.` for neither. `deceased` is 1 if the individual has died. `age` is an integer or `.`. The `kind` column can be `birth`, `planned`, `pregnancy`, `miscarriage` or `stillbirth`. The `twin` column holds `MZ:a`, `DZ:a` or `?:a`; the letter names the twin group, and sibs with the same letter are in the same group. The `adopted` column can be `in` (adopted into the family) or `out` (adopted out). `order` is the left-to-right drawing order of sibs.

A `##union` line is only needed for a partnership that cannot be derived from the file: a childless marriage, a separation, or consanguinity. Two individuals with a common child already count as partners. `##question` lines record questions that `pedkit risk --questions` answers (see section 7.5).

## 3. Validation

### 3.1. Command

```sh
pedkit validate examples/cf_family.pedx
```

The command reads the file, checks every rule and prints each problem on its own line. It does not stop at the first error; it collects all of them. For a file without problems the output is a single line:

```
6 individuals, 2 diseases · 0 errors, 0 warnings
```

The command exits with code 1 if there is at least one error and 0 otherwise, so it can be used in scripts, for example `pedkit validate family.pedx && echo valid`.

The following output was produced by deliberately corrupting three lines of `cf_family.pedx`: the sex of I-1 was set to 2, the father of II-1 was written as `I-9`, and the D1 value of II-2 was set to 5:

```
ERROR E5 (line 10) [I-1]: 'I-1' is female (sex=2) but is the father of: II-2.
ERROR E4 (line 12) [II-1]: 'I-9', given as father, is not in the file.
ERROR E8 (line 13) [II-2]: Invalid value '5' in column 'D1'; allowed: 0, 1, 2, 3.
6 individuals, 2 diseases · 3 errors, 0 warnings
```

Each line reads in this order: the level (`ERROR` or `WARNING`), the rule code, the line number of the problem in the file, the individual concerned in square brackets, and an explanation. The line number is the one your editor shows; the first line is 1.

### 3.2. Rules

An ERROR finding means the file cannot be used reliably. Mode inference and risk computation refuse a file with errors and stop after listing them. Drawing is done as far as possible.

| Code | Meaning |
|---|---|
| E1 | The first line is not `##format=PEDX/...`, or there is no column header starting with `#`. |
| E2 | A required column is missing, the first five columns are out of order, the sixth column is not the primary disease, or a column name is repeated. |
| E3 | The same individual ID is used twice in a family. The line number points to the second use. |
| E4 | A father, mother or `##union` entry refers to an individual not in the file. |
| E5 | An individual listed as a father is female, or one listed as a mother is male. The line number is the parent's own line, not the child's. |
| E6 | The pedigree has a loop: someone is their own ancestor. The finding is given on the line of the first individual of the loop in the file, and everyone in the loop is listed. |
| E7 | A disease has no phenotype column, or a `Dn`/`G:Dn` column in the header has no `##disease` definition. |
| E8 | A value is outside its allowed set: invalid sex, phenotype or genotype spelling, missing or extra cells, a haploid genotype for an autosomal disease, and similar. |
| E9 | A male has a two-allele genotype (`0/1`, `1/1`) for an X-linked disease. |

A WARNING finding does not make the file unusable, but it should be checked.

| Code | Meaning |
|---|---|
| W1 | MZ twins differ in sex or in tested genotype. |
| W2 | A tested genotype contradicts the observed phenotype under the mode given in the file (with full penetrance and no phenocopy). |
| W3 | A `##question` line refers to an individual or disease that does not exist. |
| W4 | Only one parent of an individual is given; the other is 0. |

### 3.3. Interpretation

A validation without errors shows that the file is formally correct and internally consistent. Everyone's parents are in the file, sexes agree with parental roles, values are in their allowed sets.

A validation without errors does not show that the pedigree was entered correctly in the biological sense. Marking the wrong individual as affected, swapping two sibs or entering a wrong allele frequency breaks no format rule, and validation does not catch these mistakes. They are carried into every later computation unchanged.

W2 compares each tested genotype with the phenotype the mode predicts for it. It is checked only when the file gives a Mendelian mode (not `unknown` or `MF`) and only when penetrance is 1 and phenocopy 0, because with reduced penetrance or phenocopies an unaffected carrier or an affected non-carrier is possible and is not a contradiction. The absence of W2 therefore does not mean the data fit the mode in general; that question is answered by mode inference (section 6).

W4 does not count a single-parent individual as an error. In the computations the allele from the missing parent is drawn from the population frequency. This is usually a reasonable assumption, but if you know something about the missing parent (for example that they are affected), add them to the file.

## 4. Visualisation

### 4.1. Viewing in a window

```sh
pedkit-view examples/cf_family.pedx
```

The command opens the pedigree in a window. Validation findings are printed to the terminal, and the status bar shows the numbers of individuals, errors and warnings. Hovering over the status bar shows all findings as a tooltip. If the file has a disease with `mode=MF`, the window opens with two tabs: "Pedigree" shows the drawing, and "Multifactorial risk" shows the risk report and plots described in section 8.6. `File > Save as` (Ctrl+S) saves the drawing as SVG, PDF or PNG. Ctrl+Q closes the window.

If the window does not open and a message beginning `Could not load PyQt5/QtSvg` appears, PyQt5 is either not installed for the Python running the command, or the system lacks the Qt SVG package (see section 1).

### 4.2. Drawing to a file

```sh
pedkit draw examples/cf_family.pedx -o pedigree.svg
pedkit draw examples/cf_family.pedx -o pedigree.pdf
pedkit draw examples/cf_family.pedx -o pedigree.png
```

The output format follows from the file extension. SVG needs no extra package and opens in every browser. In a browser, hovering over a symbol shows that individual's phenotype and tested genotype for each disease as a tooltip; the Qt window does not show these tooltips. PDF and PNG need PyQt5 with Qt SVG. PDF is vector output and stays sharp when enlarged, which makes it the best choice for slides and assignments. PNG is produced at twice the nominal resolution. If there are validation findings the drawing is still made and the findings are printed to the terminal.

The SVG uses only the SVG 1.2 Tiny feature set (no clip paths, patterns or markers), so it renders the same in browsers, in Qt and in vector editors.

### 4.3. Meaning of the symbols

The drawing uses the standard symbols of the National Society of Genetic Counselors (NSGC; Bennett et al. 2008, 2022). A square is a male, a circle a female, a diamond an individual of unknown sex. An affected individual's symbol is filled with the disease's pattern. If the file has more than one disease the symbol is divided vertically into equal parts, one per disease from left to right; in `cf_family.pedx` the left half shows D1 and the right half D2. A half-filled part (lower half patterned) is an intermediate phenotype. A `?` inside a part means the phenotype for that disease is unknown, and a central black dot marks a tested heterozygous carrier (an unaffected individual with genotype `0/1`).

A line from lower left to upper right through the symbol marks a deceased individual. An arrow at the lower left with the letter `P` marks the proband; an arrow without a letter marks the consultand. A small symbol with `P` inside is a planned child or a pregnancy, a small triangle is a miscarriage, and a deceased symbol with `SB` below it is a stillbirth. Square brackets mark adoption. Partners are joined by a horizontal line. A double line marks consanguinity, and two slashes across the line mark a separation. The double line needs no `##union` line: if the partners have a common ancestor in the file, the drawing finds it by itself. Two lines leaving the same point mark twins; a horizontal bar between them marks MZ (monozygotic) twins, its absence DZ twins.

Below each symbol are the individual ID, the age if given (for a deceased individual in the form `d. 67`) and the gender label if given. Roman numerals on the left mark the generations. The legend at the lower right lists the diseases, their patterns and, if given in the file, their modes. The correspondence of each element with the NSGC recommendations is listed in `docs/validation.md`.

### 4.4. What the drawing does not mean

The horizontal position of a symbol has no genetic meaning. The layout algorithm aims to avoid crossing lines and to keep sibs below their parents. Whether an individual stands left or right depends only on the `order` value in the file and on space constraints. A partner who marries into the family is always placed to the right of the family member.

The carrier dot is drawn only for an individual whose genotype was tested and found to be `0/1`. An individual inferred to be a carrier (for example the obligate carrier parent of an affected child) gets no dot, because the drawing shows observations, not inferences. For the same reason an empty (unaffected) symbol does not say that the individual's genotype is normal; it only says the individual has no signs of the disease.

The layout was designed for pedigrees of three to five generations and fewer than about 30 individuals. For individuals with several partners, or when two separate pedigrees are joined by marriage, a partner line may pass behind another symbol. The line stays behind the symbol and does not cut through it, but it can make the drawing harder to read.

## 5. PLINK .fam export

```sh
pedkit fam examples/cf_family.pedx -o family.fam
pedkit fam examples/cf_family.pedx -o family.fam --disease D2
```

The command writes the pedigree in PLINK's `.fam` format (Purcell et al. 2007; Chang et al. 2015). Every line has exactly six tab-separated columns: family, individual, father, mother, sex and phenotype. The phenotype column comes from the disease chosen with `--disease`, or from the primary disease. For D2 the output is:

```
FAM01	I-1	0	0	1	1
FAM01	I-2	0	0	2	1
FAM01	II-1	I-1	I-2	1	1
FAM01	II-2	I-1	I-2	2	1
FAM01	II-3	0	0	1	2
FAM01	III-1	II-3	II-2	0	0
```

A missing phenotype is written as `0`, not as PLINK's other customary value `-9`. Because the PLINK phenotype is binary, the intermediate phenotype (3) is exported as a case, that is 2; for a codominant disease this conversion loses information. Genotype columns, age, death and all other fields have no place in the `.fam` format. If an individual ID contains a space, PLINK could not split it, so the command reports an error.

## 6. Mode-of-inheritance inference

### 6.1. Method

This command looks at the phenotype and genotype observations of the pedigree and evaluates by which Mendelian mode the disease is transmitted. Seven modes are tried: autosomal dominant (AD), autosomal recessive (AR), autosomal codominant (ACD), X-linked dominant (XD), X-linked recessive (XR), Y-linked (Y) and maternally inherited mitochondrial (MT).

For each mode pedkit computes the likelihood, P(data | mode): the probability that, if the mode were true, a family would turn out exactly as observed. Founders, that is individuals whose parents are not in the file, receive alleles from the population frequency under Hardy–Weinberg equilibrium. The alleles are passed down the pedigree by the rules of the mode. Under an autosomal mode each parent passes one of two copies with probability 1/2 each. Under an X-linked mode a father passes his X only to his daughters. Under the Y-linked mode the Y passes only to sons. Under the mitochondrial mode a child receives the mother's mitochondria and never the father's. MZ twins share one genotype.

The observations enter as follows. A tested genotype must match exactly. Phenotypes enter probabilistically: an affected individual carrying the disease genotype is explained with probability equal to the penetrance, an affected individual without it with probability equal to the phenocopy rate. With penetrance 1 and phenocopy 0 this comparison becomes yes/no. The likelihoods of the modes are combined with Bayes' rule to give each mode's posterior probability, P(mode | data).

Two methods compute the likelihood. The default, `exact`, sums over all genotype assignments by generalised Elston–Stewart peeling (variable elimination; Elston and Stewart 1971; Cannings, Thompson and Skolnick 1978). Its cost grows with the width of the pedigree, not with the number of individuals, so it has no size limit. The second, `monte-carlo`, simulates the family millions of times (gene dropping). Each simulation gives every individual a genotype, and is compared with the observations. Simulations that explain every observation form the supporting pool; those that fail to explain at least one form the non-supporting pool. The mean weight of the simulations estimates the likelihood without bias. The Monte Carlo method gives the same likelihoods with sampling error, and in addition lists the joint genotype configurations of the family and which individual caused each conflict. These are useful for teaching and for understanding why a mode was excluded.

### 6.2. Command

```sh
pedkit infer examples/cf_family.pedx
pedkit infer examples/cf_family.pedx --method monte-carlo --seed 1
```

The exact method takes well under a second for a small family. With `--method monte-carlo` the default is one million simulations per mode, a few seconds for a six-person family. During a Monte Carlo run the progress is written to the error stream (stderr) and the report to standard output, so `> report.txt` captures the report without the progress lines.

| Option | Purpose |
|---|---|
| `--disease D2` | Analyse another disease instead of the primary one. |
| `--method` | `exact` (default), `monte-carlo`, or `both` (both reports, one after the other). |
| `-n 5000000` | Simulations per mode. Default 1 000 000. |
| `--seed 1` | Seed of the random number generator. The same seed gives the same result. Without it each run gives slightly different numbers. |
| `--modes AR,XR` | Try only the listed modes. |
| `--proposal-freq 0.25` | Turn on importance sampling (see section 6.6). |
| `--configurations 0` | Number of genotype configurations listed per pool; 0 lists all. Default 20. |
| `--summary` | Print only the diagnosis tree and the summary table. |

### 6.3. First part of the report: hierarchical diagnosis

The report starts with the diagnosis of the inheritance. For `cf_family.pedx` this part is:

```
Hierarchical diagnosis

  Inheritance  (P = 1)
    ├─ Autosomal  (P = 1 · within parent 1 · prior 0.4286)
    │  ├─ autosomal dominant  (P = 0 · within parent 0 · prior 0.1429)
    │  ├─ autosomal recessive  (P = 1 · within parent 1 · prior 0.1429)
    │  └─ autosomal codominant  (P = 0 · within parent 0 · prior 0.1429)
    ├─ Gonosomal  (P = 0 · within parent 0 · prior 0.4286)
    │  ├─ X-linked  (P = 0 · within parent 0 · prior 0.2857)
    │  │  ├─ X-linked dominant  (P = 0 · within parent 0 · prior 0.1429)
    │  │  └─ X-linked recessive  (P = 0 · within parent 0 · prior 0.1429)
    │  └─ Y-linked  (P = 0 · within parent 0 · prior 0.1429)
    └─ Mitochondrial (maternal)  (P = 0 · within parent 0 · prior 0.1429)

  Level 1: Autosomal  (P = 1, within parent 1)
  Level 2: autosomal recessive  (P = 1, within parent 1)

Diagnosis: autosomal recessive
```

Every node of the tree carries three numbers. `P` is the probability, given the data, that this branch is correct; it is the sum of the posterior probabilities of the modes below it. "Within parent" is the share of this branch given that its parent branch is correct. For example, if Gonosomal had P = 0.5 and X-linked P = 0.4, the "within parent" value of X-linked would be 0.8, read as "if the inheritance is gonosomal, it is X-linked with probability 0.8". `prior` is the probability given to the branch before any data are seen.

The "Diagnosis" line follows the most probable branch at each level: first the most probable of autosomal, gonosomal and mitochondrial, then the most probable sub-branch within it. There is no threshold in this choice. If one branch had probability 0.51 and the other 0.49, the diagnosis would still name the 0.51 branch. The diagnosis line should therefore never be read alone, only together with the probabilities of the tree. If two branches have exactly equal probability, both are named with the note "equally probable". If no mode can explain the observations, the line reads "No mode is compatible with the observations; no diagnosis." (see section 6.8).

The prior is spread equally over the seven modes, 1/7 each. This has a side effect: the Autosomal and Gonosomal branches, with three modes each, start at 3/7, while the single-mode Mitochondrial branch starts at 1/7. With strong data the difference does not show in the result. With weak data (few individuals, few affected, no tests) the mitochondrial branch starts at a disadvantage, which should be kept in mind when reading the diagnosis.

### 6.4. Second part of the report: summary table

The Monte Carlo summary for `cf_family.pedx` is:

```
Mode-of-inheritance inference · D1 (Cystic fibrosis)
q=0.02 · penetrance=1 · phenocopy=0 · mutation rate=0 · 1,000,000 simulations/mode

Mode    Supporting  Frequency   P(data|mode)        ±SE   Prior  Posterior
AR             300   3.00e-04       3.00e-04   1.73e-05   0.143          1
AD               0          0              0          0   0.143          0
...
Posterior ranking: AR (1) > AD (0) > ACD (0) > XD (0) > XR (0) > Y (0) > MT (0)
```

The header line repeats the parameters used in the computation. Every result depends on them; if you see an unexpected value here, check the file's `##disease` line.

`Supporting` is the number of simulations compatible with the observations, and `Frequency` its share of all simulations. `P(data|mode)` is the likelihood estimate. With penetrance 1, phenocopy 0 and no importance sampling it equals the frequency. With reduced penetrance the supporting simulations are weighted, so it is smaller than the frequency. `±SE` is the standard error of the likelihood estimate; the estimate may deviate from the true value by roughly two standard errors. `Prior` and `Posterior` are the probabilities of the mode before and after the data. Rows are sorted by posterior, largest first.

In the exact report the same table shows the exact likelihood, a standard error of 0, and empty `Supporting` and `Frequency` columns, since no simulation was run. For `cf_family.pedx` the exact AR likelihood is 2.88e-04; the Monte Carlo estimate 3.00e-04 ± 1.73e-05 is within one standard error of it.

In `cf_family.pedx` the likelihood of AR is about 3e-04. This number does not mean "AR has probability 0.0003". It means: if the disease is AR, a random family turns out with exactly this distribution of phenotypes and genotypes with probability about 3 in 10 000. Seeing an affected child of a rare disease is an improbable event to begin with, so likelihoods are usually very small. The size of a likelihood says nothing by itself; what matters is the ratio between modes. The posterior column makes that comparison.

### 6.5. Third part of the report: genotype pools

Each mode has its own section. In a Monte Carlo report the supporting pool lists the genotype configurations compatible with the observations, from most to least frequent:

```
── AR ──
Supporting pool: 300 simulations, 7 distinct configurations, effective sample 300
     0.34  (      102)  I-1=0/1 I-2=0/1 II-1=1/1 II-2=0/1 II-3=0/0 III-1=0/1
    0.307  (       92)  I-1=0/1 I-2=0/1 II-1=1/1 II-2=0/1 II-3=0/0 III-1=0/0
      0.3  (       90)  I-1=0/1 I-2=0/1 II-1=1/1 II-2=0/0 II-3=0/0 III-1=0/0
  Individual genotype probabilities (given the observations):
    I-1      0/1: 1
    I-2      0/1: 1
    II-1     1/1: 1
    II-2     0/0: 0.31, 0/1: 0.69
    II-3     0/0: 0.947, 0/1: 0.053
    III-1    0/0: 0.627, 0/1: 0.357, 1/1: 0.017
```

Each configuration line is a set of genotypes assigned jointly to every individual in the family. The first number is the probability of that set within the pool, the number in parentheses how many simulations produced it. The `effective sample` is the number of independent simulations the weighted simulations are worth. With equal weights it equals the number of simulations. If it is small (below a few hundred), the shares in the pool are noisy.

The individual genotype probabilities give, for the given mode and observations, each individual's genotype distribution. This table shows the consequences of a mode in the family most clearly. In `cf_family.pedx` under AR, I-1 and I-2 are `0/1` with probability 1: two healthy parents of an affected child are obligate carriers. II-2 is a carrier with probability 0.69. The exact value is 2/3, which the exact report prints; the difference comes from the sampling noise of a 300-simulation pool. These values are also the basis of the Bayes risk (section 7).

The non-supporting pool section gives three pieces of information. Its prior probability is the probability that, if the mode were true, a random family would contradict the observations. The line "Individuals causing the conflict" shows for what share of the non-supporting simulations each individual's observation could not be explained:

```
  Individuals causing the conflict (share): II-1: 1, II-3: 4.2e-04, I-2: 4.0e-04, ...
```

This line is the most useful information for understanding why a mode was excluded. In `cf_family.pedx` the share of II-1 is 1 under every mode, so the conflict almost always comes from the affected II-1. Under the X, Y and MT modes the reason is that II-1 is male and his genotype `G:D1=1/1` has two alleles. A male has a single X for an X-linked gene, and a two-allele test result is incompatible with these modes. This single test result excludes the X, Y and MT modes on its own. The shares can add up to more than 1, because several individuals can conflict in the same simulation. The section ends with the non-supporting configurations and the genotype distribution in this pool. The distribution of the non-supporting pool largely reflects the population distribution and carries less information for interpretation than the supporting pool.

### 6.6. Rare alleles and importance sampling

For a rare disease the supporting pool stays small. In `cf_family.pedx` only 300 of one million simulations supported AR, because two founders being carriers at the same time has probability (2q)², about 0.0016. The option `--proposal-freq 0.25` draws the founder alleles with frequency 0.25 instead of the true frequency and reweights each simulation by the probability it would have had under the true frequency. This is called importance sampling. The estimate stays unbiased, the supporting pool can grow a hundredfold and the standard error shrinks.

With importance sampling on, the `Supporting` and `Frequency` columns count the simulations themselves, not probabilities in the real population. The columns to read as probabilities are `P(data|mode)` and `Posterior`. Without importance sampling the two kinds of number coincide; with it they separate. Importance sampling helps most in small and medium pedigrees. In large pedigrees with many founders the weights can become very uneven and the effective sample small; there the exact method should be preferred.

### 6.7. Exact computation and enumeration

The exact method is the default because it has no sampling error and no size limit. Its results are what the Monte Carlo estimates converge to, so comparing the two is a check on whether the simulation was long enough. A Monte Carlo likelihood of 0 does not mean a mode is excluded; it means that no compatible simulation was seen. In `mitochondrial.pedx` the exact computation gives AR a likelihood of about 2.4e-09, while one million simulations give 0: the expected number of compatible simulations is 0.002. An exact likelihood of 0, in contrast, really means impossible.

A third, independent exact method is available from Python and for the Bayes risk: enumeration sums over every joint genotype assignment one by one (`pedkit.inference.enumeration_likelihood`, `pedkit risk --method enumeration`). Its cost grows exponentially with the number of individuals, so it is limited to 14 individuals. It shares no code with peeling except the genetic model, which makes it a useful cross-check; the test suite compares the two on every example (see `docs/validation.md`).

### 6.8. What mode inference does and does not say

A posterior of 1 does not say that the mode is certainly correct. It says: among the seven Mendelian modes tried, with the parameters in the file (q, penetrance, phenocopy, mutation rate), only this mode can explain the observations. If the true mode is not one of the seven, the tool cannot know it. Imprinting, anticipation, digenic inheritance, age-dependent penetrance and mitochondrial heteroplasmy are not modelled.

The parameters change the result. In `cf_family.pedx` AD comes out as 0, because an affected child of two healthy parents is impossible under AD with penetrance 1 and no new mutation. If the same file has penetrance 0.8 or a mutation rate of 0.001, AD is no longer 0. When you see a mode excluded, ask which parameter assumption the exclusion rests on.

A single test result can decide the whole inference. If the `G:D1` cell of II-1 in `cf_family.pedx` is set to `.` (the test result removed), the result changes completely: XR becomes about 0.985 and AR about 0.015. The reason is that a male inheriting an X-linked disease from a carrier mother (about 2q) is much more likely than inheriting a recessive disease from two carrier parents (about q²). This difference between the tested and untested results is not a fault of the tool; it shows the value of the test. To see why a mode was chosen, look at the "Individuals causing the conflict" line and at the tested individuals.

A high posterior from one small family is not population-level evidence. Transmission in a single family can be compatible with several modes. A clinical or publication-level decision on a mode needs several families and molecular data evaluated together.

## 7. Bayes risk

### 7.1. Method

The Bayes risk is the probability that a particular individual is in a particular state: being a carrier, being affected, or having a particular genotype. The computation builds the classical table used in genetic counselling. Its rows are the possible genotypes of the individual, and four numbers are computed for each row.

The prior is the probability of the genotype given the individual's own observations and the information from ancestors and collateral relatives. The conditional is the probability of the observations of the individual's descendants and of their married-in partners, given the genotype. The joint is the product of prior and conditional. The posterior is the joint divided by the sum over rows; it is the probability of the genotype given all data. The split of the data into prior and conditional information is automatic, and the start of the report lists which individuals give the conditional information.

If the disease's mode is given in the file (such as `mode=XR`) the computation uses only that mode. With `mode=unknown` a separate table is built for each mode and the risk is averaged with the posterior probabilities of the modes, which are the same as in mode inference.

### 7.2. Command

```sh
pedkit risk examples/dmd_bayes.pedx --individual II-2 --event carrier
```

| Option | Purpose |
|---|---|
| `--individual II-2` | The individual whose risk is computed. In a file with several families, write `FAM04/II-2`. |
| `--event` | The state to compute. `carrier` is being heterozygous (`0/1`); `mutant` is carrying at least one disease allele (`0/1`, `1/1` or hemizygous `1`); `affected` is being affected; `intermediate` is showing the intermediate phenotype (codominant); `genotype:1/1` is having a particular genotype. Default `carrier`. |
| `--disease D2` | Another disease instead of the primary one. |
| `--modes XR` | Compute only under this mode, overriding the `mode` field of the file. Several modes can be given with commas (`AR,XR`). |
| `--method` | `exact` (peeling, no size limit), `enumeration` (at most 14 individuals; an independent check), `monte-carlo`, or `auto` (= exact). |
| `-n`, `--seed`, `--proposal-freq` | Monte Carlo settings as in section 6.2. |
| `--questions` | Answer every `##question` line of the file (section 7.5). |

### 7.3. Output and interpretation

`examples/dmd_bayes.pedx` is a textbook example for Duchenne muscular dystrophy. II-2's brother II-1 is affected, and II-2 has two healthy sons. The question is the probability that II-2 is a carrier.

```
Bayes risk · II-2 · D1 (Duchenne muscular dystrophy)
Event: heterozygous carrier (0/1) · method: exact
Prior information: the individual's own observations, ancestors and collateral relatives
Conditional information: II-3, III-1, III-2, III-3

── XR · X-linked recessive ──
  Genotype                Prior         Conditional               Joint           Posterior  P(event|g)
  0/0                 1/2 (0.5)              0.9999                 0.5           4/5 (0.8)           0
  0/1                 1/2 (0.5)                0.25               0.125           1/5 (0.2)           1
  P(heterozygous carrier (0/1) | data, XR) = 1/5 (0.2)

Bayes risk: P(heterozygous carrier (0/1)) = 1/5 (0.2)
```

The prior is 1/2. II-1 being affected makes their mother I-2 an obligate carrier, and II-2 received one of her mother's two X chromosomes with probability 1/2 each. The conditional column comes from the two healthy sons. If II-2 is a carrier each son is healthy with probability 1/2, both with probability 1/4. If she is not a carrier the sons are healthy in any case. The joint column is 1/2 × 1/4 = 1/8 and 1/2 × 1 = 1/2. The posterior divides them by their sum: (1/8) / (1/8 + 1/2) = 1/5. The two healthy sons lowered II-2's carrier probability from 1/2 to 1/5.

The value 0.9999 in the 0/0 row comes from the partner II-3 being part of the conditional information. II-3 is healthy with probability 1 − q, which does not depend on II-2's genotype, so it is the same in both rows and cancels in the posterior. Textbooks leave this term out because it has no effect on the result. Showing it in the table has one benefit: it shows who enters the conditional information. The conditional of the 0/1 row is also really 0.25 × 0.9999, but rounds to 0.25 at four digits.

The `P(event|g)` column is the probability of the asked event if the row's genotype is true. For the carrier question it is 1 in the 0/1 row and 0 elsewhere. For the affected question it is the penetrance in rows with the disease genotype and the phenocopy rate elsewhere. The risk is the sum of the products of the posterior and this column.

When a value corresponds to a simple fraction, the fraction is printed with the decimal in parentheses. Fractions are printed only for exact methods. Monte Carlo results are printed only as decimals, because an estimate such as 0.69 has no exact ratio 69/100.

### 7.4. Mode averaging

If no mode is given, the report first lists the mode weights:

```
Mode weights P(mode | data):
  AR                        1   P(event | data, AR) = 2/3 (0.6667)
  AD                        0   P(event | data, AD) = —
  ...
```

Each line shows a mode's posterior probability and the individual's risk if that mode is true. The final Bayes risk is the sum of their products. A `—` means that the mode is incompatible with the observations and the risk under it is undefined. Tables of modes with weight 0 are not printed.

The averaged risk carries the uncertainty about the mode into the result. If II-1's test result is removed from `cf_family.pedx`, II-2's carrier risk is 1/2 under XR and 2/3 under AR, and the averaged risk is their average weighted by about 0.985 and 0.015. If the mode is known for certain, write it in the file or give it with `--modes`; otherwise the tool gives every mode its share.

### 7.5. Answering the questions in a file

```sh
pedkit risk examples/cf_family.pedx --questions
```

The command answers each `##question` line in turn. A `type=mode` question produces a mode inference report (with the hierarchical diagnosis), a `type=risk` question a Bayes risk report. Answers are separated by rules. For `cf_family.pedx` the answer to Q1 is autosomal recessive, to Q2 (II-2 is a carrier) 2/3, and to Q3 (III-1 is affected) 1/153. The last value equals 2/3 (II-2 a carrier) × 2/51 (II-3 a carrier, that is 2q/(1+q)) × 1/4 (an affected child of two carriers). If a question cannot be answered (for example because it refers to a missing individual) the reason is printed and the other questions are still answered.

### 7.6. What the Bayes risk does and does not say

The risk is the probability that the individual is in the state, under the given mode and parameters, with every observation of the file taken into account. If the mode, allele frequency or penetrance is entered wrongly, the risk is wrong and the tool cannot notice. In the DMD example the mutation rate was entered as 0, so I-2 was taken as an obligate carrier. In real DMD about a third of isolated cases are new mutations. With a realistic `mutation_rate` I-2 is no longer certainly a carrier and II-2's risk drops below 1/5.

For an individual whose phenotype was observed, the "affected" risk is the observation itself: 1 for II-1, 0 for II-2. The tool does not compute the probability that someone healthy today develops a late-onset disease later, because age-dependent penetrance is not modelled. For such a question use the carrier or `mutant` event and interpret the result separately with the penetrance.

The carrier risk is the probability of heterozygosity. For an X-linked disease a male is hemizygous and cannot be heterozygous, so the `carrier` probability of a male is always 0. To ask whether a male carries the disease allele, use `--event mutant` or `--event genotype:1`.

For a planned child or a pregnancy of unknown sex, the X-linked risk is computed with each sex having probability 1/2. In the DMD example the 1/20 for III-3 is the average of "1/10 if a boy" and "0 if a girl". If the sex is known, write it in the file; the result changes.

A risk computed by Monte Carlo is an estimate. The report states the number of compatible simulations in each mode's heading, and the risk should be judged with that number in mind. A value of 0.69 from a few hundred compatible simulations is an ordinary deviation when the true value is 0.667. The exact method has no such error and no size limit, and is the default.

## 8. Multifactorial risk

### 8.1. Model

Diseases such as cleft lip, pyloric stenosis, neural tube defects or congenital hip dislocation do not arise from a single gene but from the sum of many small-effect genes and the environment. A Mendelian genotype is not defined for them, so the Bayes risk of section 7 cannot be used. Instead pedkit uses Falconer's liability-threshold model (Falconer 1965). Every individual has a normally distributed liability with population mean 0 and standard deviation 1, and an individual whose liability exceeds a threshold is affected. The threshold follows from the prevalence: for a prevalence K it is the value T that cuts off the upper K tail of the standard normal distribution (T = 2.326 for K = 0.01). A share h² of the liability is heritable and is shared between relatives in proportion to their relationship; the liability correlation of sibs and of parent–child pairs is h²/2, of second-degree relatives h²/4, and of MZ twins h². The remaining 1 − h² is the individual's own environment.

An individual's risk is the probability that their liability exceeds the threshold, given which relatives are affected or unaffected. An affected relative shows that the liability shared by the family is above average and shifts the individual's liability distribution toward the threshold. Unaffected relatives pull it back a little. The computation uses the information of the whole family jointly.

pedkit computes the same risk by two independent methods and reports them side by side. The first is Elston–Stewart peeling. A continuous liability cannot be peeled directly, so the heritable part is represented by the hypergeometric polygenic model with n loci (Cannings, Thompson and Skolnick 1978; Fernando, Stricker and Elston 1994): each parent passes n of their 2n alleles to a child without replacement, which preserves the parent–child resemblance and the Mendelian sampling variance. Peeling uses generalised variable elimination rather than nuclear families, so consanguinity loops are handled correctly. At finite n the tail of the distribution is thinner than the true normal tail, and for a rare disease the risk is decided in exactly that tail. The risk is therefore computed with n and with 2n loci, and Richardson extrapolation gives the value for infinitely many loci. The second method is MCMC. The liabilities of all individuals are drawn by Gibbs sampling from the multivariate normal distribution built from the kinship matrix, truncated to the intervals the observed phenotypes allow (Albert and Chib 1993). The target's risk is the average, over the steps, of the probability that the target exceeds the threshold given the other individuals' liabilities (a Rao–Blackwell estimate). Sixty-four independent chains run at once; the standard error comes from the spread of the chain means, and convergence is measured by the Gelman–Rubin R̂ (Gelman and Rubin 1992).

### 8.2. Defining the disease in the file

The disease is marked with `mode=MF` on its `##disease` line, and the two parameters of the model go on the same line:

```
##disease=<ID=D1,name="Cleft lip",mode=MF,prevalence=0.001,heritability=0.76,pattern=solid>
```

For diseases whose prevalence differs between the sexes, write `prevalence_male` and `prevalence_female` instead of `prevalence`. Each sex then has its own threshold. If these parameters are missing from the file, or you want to try other values, give them on the command line with `--prevalence`, `--heritability`, `--prevalence-male` and `--prevalence-female`; they override the file. The intermediate phenotype (3) is not defined in this model; the phenotype column may only hold 0, 1 and 2.

### 8.3. Command

```sh
pedkit mf examples/multifactorial.pedx --individual III-2
pedkit mf examples/multifactorial.pedx --individual III-2 --plots plots/
pedkit mf examples/multifactorial.pedx --questions --plots plots/ --format svg
pedkit mf examples/multifactorial.pedx --individual III-2 --heritability 0.5 --method es
```

`--method` can be `both` (default), `es` or `mcmc`. `--loci` is the n of Elston–Stewart (default 20; the computation uses 20 and 40 loci). `--chains`, `--samples` and `--burn-in` set the number of MCMC chains, the steps kept per chain and the discarded burn-in steps. `--seed` makes the MCMC reproducible. `--questions` answers every `type=risk` question line of the file. For a file whose disease has `mode=MF`, `pedkit risk --questions` also answers the risk questions with this model, but `pedkit risk --individual` called directly gives an error that points to this command.

### 8.4. Output and interpretation

The example file uses infantile hypertrophic pyloric stenosis. The disease is five times more common in males (prevalence 0.005 in males, 0.001 in females), and its heritability is taken as 0.8. The mother (II-2) is affected, III-2 is a planned son and III-3 a planned daughter. The first part of the report for III-2 is (the MCMC values depend on random sampling and change in the last digits between runs; `--seed` fixes them):

```
Multifactorial risk · III-2 · D1 (Infantile hypertrophic pyloric stenosis)
Model: liability-threshold (Falconer) · h² = 0.8
  male: K = 0.005, threshold T = 2.576 · female: K = 0.001, T = 3.090
Individual: III-2 · sex male · phenotype unknown

Risk of being affected given the relatives' phenotypes:
  Elston–Stewart (HPM, Richardson n→∞): 0.0816
      n=20 loci: 0.0737 · n=40 loci: 0.0777
  MCMC (Gibbs, 64 chains × 4000 steps, burn-in 1000): 0.0818
      standard error 2.384e-04 · count estimate (L > T) 0.0823
      R̂ liability 1.0001 · R̂ risk 1.0002
  Difference of the two methods: +0.8 standard errors
  Population risk (prevalence): 0.005
  Relative risk (λ = risk / prevalence): 16.33
```

The son's risk is about 0.082, sixteen times that of males in the population. For a daughter of the same mother the risk is 0.025. That is 25 times her population risk, but lower in absolute terms than her brother's, because the threshold of girls is higher. The Carter effect follows from the same model: an affected individual of the less often affected sex must carry a higher liability because their threshold is higher, so their relatives' risk is higher than that of the relatives of an affected male.

The model lines at the start give the parameters used: heritability, the prevalence of each sex and the threshold derived from it. When values from the command line or the viewer override the file, the values shown here are the ones used. HPM stands for hypergeometric polygenic model. The n=20 and n=40 lines are the inputs of the Richardson extrapolation; the larger their difference, the larger the finite-locus error, and the extrapolated value removes its first-order part. The first line under MCMC gives the standard error of the estimate. The "Difference of the two methods" line gives how many standard errors the Elston–Stewart result lies from the MCMC estimate. R̂ is the ratio of between-chain to within-chain variability and is close to 1 when the chains sample the same distribution. The count estimate is the share of sampled liabilities above the threshold; it estimates the same quantity with more noise and is shown only for comparison.

The second part of the report is the table of individuals:

```
  Indiv.  Sex    Phenotype         φ   ES risk*  L mean   L 90% interval
  I-1     male   unaffected    0.125      0.081   +1.16   [-0.23, +2.35]
  II-2    female affected       0.25  8.026e-04   +3.35   [+3.10, +3.85]
  III-2   male   unknown         0.5     0.0816   +1.30   [-0.21, +2.80]
```

φ is each individual's kinship coefficient with the target (0.5 for the individual itself and for an MZ co-twin, 0.25 for a sib or parent, 0.125 for a second-degree relative). An individual of unknown sex has `?` in the Sex column. The "ES risk" column is each individual's risk from the relatives if their own phenotype were unknown; for the affected mother, for example, it is her probability of being affected given that her parents and brother are unaffected. "L mean" is the MCMC posterior mean of the liability given all data, and "L 90% interval" the 5 % and 95 % quantiles of its posterior distribution. The lower end of an affected individual's interval equals their threshold. With Elston–Stewart only (`--method es`), the last two columns show `—`. Report lines do not exceed 76 characters, so they fit an 80-column terminal and the viewer's report panel without scrolling.

The target's own phenotype is ignored if it was observed. The risk computed for an affected individual is not 1 but the value that would be expected if their status were unknown. For a target of unknown sex the thresholds of both sexes are averaged with weight 1/2 each.

### 8.5. Plots

`--plots FOLDER` draws five plots. `liability` overlays the target's liability distribution on the population N(0, 1) distribution; the shaded area right of the threshold is the risk. The Elston–Stewart curve and the MCMC histogram are two estimates of the same distribution. `heritability` shows how the risk changes with h² at fixed prevalence; at h² = 0 the risk falls to the population risk. `prevalence` shows on a log scale how the risk changes with the prevalence at fixed h²; the dashed line is risk = prevalence, and the vertical distance to the curve is the logarithm of the relative risk. `relatives` shows the "ES risk" column of the table together with the prevalence in each individual's own sex. `mcmc` has two panels: above, the cumulative mean of the risk estimate with a ±2 standard error band and the Elston–Stewart value; below, the traces of the target's liability in the first three chains over the first 400 steps, with the threshold as a dashed line. If the chains mix and move back and forth across the threshold, the sampler explores the distribution; a chain staying at a separate level for a long time also pushes R̂ away from 1. `--format` can be png (default), svg or pdf. The plots need matplotlib.

### 8.6. Multifactorial risk in the viewer

```sh
pedkit-view examples/multifactorial.pedx
```

If the file has a disease with `mode=MF`, the viewer window shows a "Multifactorial risk" tab next to the pedigree. The top row of the tab selects the disease, the individual, the method and h². The "Compute" button shows the report on the left and the five plots in separate sub-tabs on the right; the line at the bottom of the window summarises the risk, the population risk and λ. The individual box opens with the individual of the file's risk question, or else with the first individual of unknown phenotype. The h² box starts with the value from the file. Changing it and computing again shows one point of the heritability curve with its full report. If the file gives no heritability, the box starts at 0.5 and a note at the bottom says so. The prevalence cannot be changed from the window; to try another prevalence use the `--prevalence` options of the command line.

When a computation finishes, the report panel is widened so that its longest line fits without horizontal scrolling, and the rest of the space goes to the plots. When the window is enlarged only the plot area grows. The divider between the two panels can be dragged to change the split by hand. The matplotlib toolbar above each plot has zoom, pan and save buttons for that plot; the "Save plots…" button writes all five to a chosen folder as png, svg or pdf.

The computation runs in a separate thread, so the window does not freeze; for the example family both methods with plots take about 4 seconds, Elston–Stewart alone about 1 second. "Compute" cannot be pressed again until the computation has finished. If it ends with an error (for example a missing parameter), the message appears in the bottom line and the previous results stay in place. If numpy or matplotlib is missing, the tab is not shown, the reason is printed to the terminal, and the pedigree is still displayed.

### 8.7. What the multifactorial risk does and does not say

The computation assumes that liability is fully additive, that relatives share no common environment, that liability is distributed the same way in both sexes, and that only the threshold differs by sex. If the common environment of sibs raised in the same home is absorbed into the heritability, h² is entered too high and the risk is inflated. The single-threshold model does not distinguish severity; it puts severe and mild cases in the same "affected" category. Founders are assumed to be unrelated and drawn at random from the population; the selection bias of having found the family through an affected individual (ascertainment) is not corrected.

Mode inference (section 6) does not test the multifactorial mode. It only compares the seven Mendelian modes with each other, and for a file with `mode=MF` it still reports one of them as most probable; in the example family this is autosomal recessive, because a single affected daughter of healthy parents is most easily explained that way. This result does not show that the disease is Mendelian. That a disease is multifactorial is learned not from one family but from twin concordance, from recurrence risks that fall quickly with the degree of relationship, and from sex differences; the tool takes this knowledge from the user as `mode=MF`.

The result is sensitive to the parameters. In the example family, lowering h² from 0.8 to 0.6 lowers the son's risk from about 0.08 to about 0.05; the heritability plot shows this sensitivity directly. Published h² estimates for the same disease vary between studies, so the reported risk holds if the entered parameters are correct. Empirical recurrence risks in textbooks (for example about 4 % for cleft lip after one affected sib) were counted in observed families, and comparing them with the model result is a way to check whether the parameters are reasonable.

## 9. Use from Python

Everything the command line does can be done from Python. This is useful for using results in other computations or processing many files in a loop.

```python
from pedkit import read, validate, write, write_fam
from pedkit.vis import draw, save
from pedkit.inference import infer_mode, infer_mode_exact
from pedkit.risk import bayes_risk, answer_questions
from pedkit.mf import mf_risk

ped = read("examples/cf_family.pedx")       # a Pedigree object
for finding in validate(ped):               # list of Finding; empty if no problems
    print(finding.level, finding.rule, finding.line, finding.message)

save(ped, "cf_family.pdf")                  # .svg, .pdf or .png
svg_text = draw(ped)                        # the SVG as a string

result = infer_mode_exact(ped, "D1")
print(result.report())
print(result.diagnosis.diagnosis)           # "autosomal recessive"
ar = result.results[result.most_likely.mode]
print(ar.posterior, ar.likelihood, ar.supporting.marginal["II-2"])

mc = infer_mode(ped, "D1", n=1_000_000, seed=1)   # Monte Carlo with genotype pools

risk = bayes_risk(ped, "II-2", "carrier")
print(risk.risk)                            # 0.666...
for m in risk.modes:
    for row in m.rows:
        print(m.mode.value, row.genotype, row.prior, row.conditional, row.posterior)

mf = mf_risk(read("examples/multifactorial.pedx"), "III-2", seed=1)
print(mf.es.risk, mf.mcmc.risk, mf.mcmc.standard_error, mf.relative_risk)
mf.plots("plots/", "svg")                   # returns the paths written
print(mf.report())
```

The plots can also be obtained as matplotlib `Figure` objects without writing files. These objects are not tied to pyplot; they can be shown in a Jupyter notebook or embedded in another Qt window with `FigureCanvasQTAgg`:

```python
from pedkit.mf.plots import figures
fig = figures(mf)["heritability"]           # liability, heritability, prevalence, relatives, mcmc
fig.savefig("heritability.pdf")
```

A file can be changed in Python and saved. Unchanged lines are written back byte for byte; only the changed lines are regenerated:

```python
ind = ped.individual("II-2")
ind.age = 28
write(ped, "cf_family_updated.pedx")
write_fam(ped, "cf_family.fam", disease="D1")
```

## 10. Common error messages

`ModuleNotFoundError: No module named 'numpy'` (or `'PyQt5'`, `'matplotlib'`) means the optional dependency is missing for the Python running the command. Install the matching extra (section 1), or activate the virtual environment it was installed in.

`The pedigree has errors; inference is not possible` appears in mode inference and risk computation. The validation errors are listed below the message. Fix the file first with `pedkit validate`.

`'IV-9' is not in the pedigree.` means the ID given with `--individual` is not in the file. IDs are case-sensitive: `ii-2` and `II-2` differ.

`Enumeration is limited to 14 individuals` means `--method enumeration` was used for a large family. Use `--method exact`, which has no limit.

"No mode is compatible with the observations; no diagnosis." means none of the seven modes can explain the observations with the parameters in the file. This is not a program error. Either the data contain an entry mistake (affected status on the wrong individual, a wrong test result), or the disease has reduced penetrance, phenocopies or new mutations. Check the data first, then the parameters of the `##disease` line. The message also appears for a file with intermediate phenotypes (3) when `--modes` does not include ACD, because only the codominant mode can explain an intermediate phenotype.

In the multifactorial risk the most common messages concern the parameters. `No heritability for D1` or `No prevalence for D1` means the `##disease` line lacks `heritability=` or `prevalence=` (or `prevalence_male=` with `prevalence_female=`); add the field to the file or give it on the command line. `Heritability must be in [0, 1)` appears for h² = 1, because the environmental variance cannot be zero in the model. `The intermediate phenotype (3) is not defined in the liability-threshold model` means the phenotype column contains a 3. `D1 is multifactorial (mode=MF); a Bayes risk over Mendelian genotypes is not defined` means `pedkit risk --individual` was called for an MF disease; the same computation is done by `pedkit mf`. Conversely, `a multifactorial risk needs mode=MF` means `pedkit mf` was called for a disease with a Mendelian mode. `The observations are impossible under these parameters (likelihood 0)` means that a very small prevalence or a very low heritability together with many affected relatives gave a numerically zero probability; review the parameters.

## 11. Example files

`examples/cf_family.pedx` is a small three-generation family with two diseases. It contains all basic elements of the format and three `##question` lines. `examples/codominant.pedx` shows a sickle-cell-like codominant disease with intermediate phenotype (3) values and gives an autosomal codominant diagnosis. `examples/mitochondrial.pedx` shows a disease that affected mothers pass to all their children and affected fathers to none, and gives a mitochondrial diagnosis. `examples/dmd_bayes.pedx` is the DMD Bayes example of section 7.3 with two `##question` lines. `examples/multifactorial.pedx` is the pyloric stenosis example of section 8.4; it defines the disease with `mode=MF` and has two risk questions, for a son and a daughter. Copying one of these files and editing it leads to fewer mistakes than writing a new one from scratch.

## 12. Assumptions and limits of the model

The computations of sections 6 and 7 assume a single gene with two alleles, one normal and one disease allele (the assumptions of the multifactorial model are in section 8.7). Founder genotypes are drawn from a population in Hardy–Weinberg equilibrium, so for founders from closed, consanguineous communities the carrier probability is underestimated. Penetrance does not depend on age or sex. In the X-linked dominant mode X inactivation in heterozygous females is not taken into account. In the mitochondrial mode there is no heteroplasmy: an individual either carries the mutant mitochondria or does not. New mutation acts only from the normal to the disease allele. MZ twins share one genotype (and, in the multifactorial model, one polygenic value); DZ twins and twins of unknown zygosity are treated as ordinary sibs.

These assumptions are the same as those of the pedigree problems of medical genetics courses, and the tool was designed at this scale. A risk intended for real clinical counselling must not be used before each of these assumptions has been checked for the disease concerned.

## References

Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and polychotomous response data. *Journal of the American Statistical Association* 88, 669–679.

Bennett, R. L., French, K. S., Resta, R. G. and Doyle, D. L. (2008). Standardized human pedigree nomenclature: update and assessment of the recommendations of the National Society of Genetic Counselors. *Journal of Genetic Counseling* 17, 424–433.

Bennett, R. L., French, K. S., Resta, R. G. and Austin, J. (2022). Practice resource-focused revision: Standardized pedigree nomenclature update centered on sex and gender inclusivity: A practice resource of the National Society of Genetic Counselors. *Journal of Genetic Counseling* 31, 1238–1248.

Cannings, C., Thompson, E. A. and Skolnick, M. H. (1978). Probability functions on complex pedigrees. *Advances in Applied Probability* 10, 26–61.

Elston, R. C. and Stewart, J. (1971). A general model for the genetic analysis of pedigree data. *Human Heredity* 21, 523–542.

Falconer, D. S. (1965). The inheritance of liability to certain diseases, estimated from the incidence among relatives. *Annals of Human Genetics* 29, 51–76.

Fernando, R. L., Stricker, C. and Elston, R. C. (1994). The finite polygenic mixed model: an alternative formulation for the mixed model of inheritance. *Theoretical and Applied Genetics* 88, 573–580.

Gelman, A. and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science* 7, 457–472.
