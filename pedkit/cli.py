"""Command line interface.

    pedkit validate FILE.pedx
    pedkit draw FILE.pedx -o pedigree.svg            (.svg; .pdf/.png need PyQt5 + QtSvg)
    pedkit fam FILE.pedx -o family.fam [--disease D2]
    pedkit infer FILE.pedx [--method exact|monte-carlo]
    pedkit risk FILE.pedx --individual II-2 --event carrier
    pedkit risk FILE.pedx --questions
    pedkit mf FILE.pedx --individual III-2 [--plots plots/]
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .io import read, write_fam
from .model import Level, Mode
from .validate import validate

SEPARATOR = "\n\n" + "═" * 72 + "\n\n"


def _modes(text: str | None) -> list[Mode] | None:
    if not text:
        return None
    try:
        return [Mode(m.strip()) for m in text.split(",")]
    except ValueError as e:
        raise SystemExit(f"Unknown mode in {text!r}; valid: "
                         + ", ".join(m.value for m in Mode)) from e


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="pedkit",
        description="Pedigree validation, NSGC drawing, mode-of-inheritance inference, "
                    "Bayes risk and multifactorial risk for PEDX files.")
    ap.add_argument("--version", action="version", version=f"pedkit {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate", help="validate a file and list the findings")
    p.add_argument("file")

    p = sub.add_parser("draw", help="draw the pedigree (.svg; .pdf/.png need PyQt5 + QtSvg)")
    p.add_argument("file")
    p.add_argument("-o", "--output", required=True, help="output file (.svg, .pdf, .png)")

    p = sub.add_parser("fam", help="write a PLINK .fam file")
    p.add_argument("file")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--disease", help="phenotype column (default: primary disease)")

    p = sub.add_parser("infer", help="mode-of-inheritance inference")
    p.add_argument("file")
    p.add_argument("--disease", help="disease ID (default: primary disease)")
    p.add_argument("--method", choices=["exact", "monte-carlo", "both"], default="exact",
                   help="exact peeling, Monte Carlo gene dropping, or both (default: exact)")
    p.add_argument("--modes", help="comma-separated modes, e.g. AD,AR,XR (default: all Mendelian)")
    p.add_argument("-n", type=int, default=1_000_000, help="simulations per mode (1000000)")
    p.add_argument("--batch", type=int, default=100_000, help="batch size (100000)")
    p.add_argument("--seed", type=int, help="random seed (for reproducibility)")
    p.add_argument("--proposal-freq", type=float,
                   help="founder allele frequency for importance sampling (e.g. 0.25)")
    p.add_argument("--configurations", type=int, default=20,
                   help="configurations listed per pool; 0 = all (20)")
    p.add_argument("--summary", action="store_true", help="print only the summary table")

    p = sub.add_parser("risk", help="Bayes risk for an individual")
    p.add_argument("file")
    p.add_argument("--individual", help="individual whose risk is computed (e.g. II-2)")
    p.add_argument("--event", default="carrier",
                   help="carrier, mutant, affected, intermediate or genotype:0/1 (default: carrier)")
    p.add_argument("--disease", help="disease ID (default: primary disease)")
    p.add_argument("--modes", help="only these modes, comma-separated (e.g. XR or AR,XR); "
                                   "default: the mode in the file, or all modes")
    p.add_argument("--method", choices=["auto", "exact", "enumeration", "monte-carlo"],
                   default="auto")
    p.add_argument("-n", type=int, default=1_000_000, help="Monte Carlo simulations")
    p.add_argument("--seed", type=int)
    p.add_argument("--proposal-freq", type=float, help="importance-sampling frequency")
    p.add_argument("--questions", action="store_true",
                   help="answer the ##question lines of the file")

    p = sub.add_parser("mf", help="multifactorial (liability-threshold) risk")
    p.add_argument("file")
    p.add_argument("--individual", help="individual whose risk is computed (e.g. III-2)")
    p.add_argument("--disease", help="disease ID (default: primary disease)")
    p.add_argument("--prevalence", type=float, help="population prevalence K")
    p.add_argument("--heritability", type=float, help="heritability of liability h²")
    p.add_argument("--prevalence-male", type=float)
    p.add_argument("--prevalence-female", type=float)
    p.add_argument("--method", choices=["both", "es", "mcmc"], default="both")
    p.add_argument("--loci", type=int, default=20,
                   help="HPM loci n; n and 2n are computed (default 20)")
    p.add_argument("--chains", type=int, default=64)
    p.add_argument("--samples", type=int, default=4000, help="steps per chain after burn-in")
    p.add_argument("--burn-in", type=int, default=1000)
    p.add_argument("--seed", type=int)
    p.add_argument("--plots", metavar="FOLDER", help="write the risk plots to this folder")
    p.add_argument("--format", choices=["png", "svg", "pdf"], default="png")
    p.add_argument("--questions", action="store_true",
                   help="answer the type=risk ##question lines of the file")
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = _parser()
    a = ap.parse_args(argv)
    try:
        ped = read(a.file)
    except OSError as e:
        print(e, file=sys.stderr)
        return 1
    findings = validate(ped)
    n_err = sum(f.level is Level.ERROR for f in findings)

    if a.command == "validate":
        for f in findings:
            print(f)
        print(f"{len(ped.individuals)} individuals, {len(ped.diseases)} diseases · "
              f"{n_err} errors, {len(findings) - n_err} warnings")
        return 1 if n_err else 0

    for f in findings:
        print(f, file=sys.stderr)
    try:
        return _COMMANDS[a.command](ped, a, ap)
    except (KeyError, ValueError, RuntimeError) as e:
        print(e, file=sys.stderr)
        return 1


def _draw(ped, a, ap) -> int:
    from .vis import save
    save(ped, a.output)
    print(f"Written: {a.output}", file=sys.stderr)
    return 0


def _fam(ped, a, ap) -> int:
    write_fam(ped, a.output, a.disease)
    print(f"Written: {a.output}", file=sys.stderr)
    return 0


def _infer(ped, a, ap) -> int:
    from .inference import MENDELIAN_MODES, infer_mode, infer_mode_exact
    modes = _modes(a.modes) or list(MENDELIAN_MODES)
    out = []
    if a.method in ("exact", "both"):
        out.append(infer_mode_exact(ped, a.disease, modes).report(details=not a.summary))
    if a.method in ("monte-carlo", "both"):
        def progress(mode: Mode, done: int, total: int) -> None:
            print(f"\r{mode.value}: {done:,}/{total:,}", end="", file=sys.stderr, flush=True)
            if done == total:
                print(file=sys.stderr)
        r = infer_mode(ped, a.disease, modes, n=a.n, batch=a.batch, seed=a.seed,
                       proposal_freq=a.proposal_freq, progress=progress,
                       max_configurations=a.configurations or None)
        out.append(r.report(details=not a.summary))
    print(SEPARATOR.join(out))
    return 0


def _risk(ped, a, ap) -> int:
    from .risk import answer_questions, bayes_risk
    if a.questions:
        answers = answer_questions(ped, n=a.n, seed=a.seed, method=a.method,
                                   proposal_freq=a.proposal_freq)
        if not answers:
            print("The file has no ##question lines.", file=sys.stderr)
            return 1
        print(SEPARATOR.join(x.report() for x in answers))
        return 1 if any(x.error for x in answers) else 0
    if not a.individual:
        ap.error("risk needs --individual or --questions")
    r = bayes_risk(ped, a.individual, a.event, a.disease, _modes(a.modes), method=a.method,
                   n=a.n, seed=a.seed, proposal_freq=a.proposal_freq)
    print(r.report())
    return 0


def _mf(ped, a, ap) -> int:
    from .mf import mf_risk
    common = dict(method=a.method, prevalence=a.prevalence, heritability=a.heritability,
                  prevalence_male=a.prevalence_male, prevalence_female=a.prevalence_female,
                  loci=a.loci, chains=a.chains, samples=a.samples, burn_in=a.burn_in,
                  seed=a.seed)
    if a.questions:
        jobs = [(q.individual, q.disease, q.id) for q in ped.questions
                if q.type == "risk" and q.individual]
        if not jobs:
            print("The file has no type=risk ##question line with an individual.",
                  file=sys.stderr)
            return 1
    elif a.individual:
        jobs = [(a.individual, a.disease, None)]
    else:
        ap.error("mf needs --individual or --questions")

    reports, status = [], 0
    for individual, disease, qid in jobs:
        title = f"Question {qid}\n\n" if qid else ""
        try:
            r = mf_risk(ped, individual, disease, **common)
        except (KeyError, ValueError) as e:
            reports.append(f"{title}Could not be answered: {e}")
            status = 1
            continue
        text = r.report()
        if a.plots:
            folder = a.plots if len(jobs) == 1 else f"{a.plots}/{individual.replace('/', '_')}"
            try:
                paths = r.plots(folder, a.format)
                text += "\n\nPlots:\n" + "\n".join(f"  {p}" for p in paths)
            except RuntimeError as e:
                text += f"\n\nCould not draw the plots: {e}"
                status = 1
        reports.append(title + text)
    print(SEPARATOR.join(reports))
    return status


_COMMANDS = {"draw": _draw, "fam": _fam, "infer": _infer, "risk": _risk, "mf": _mf}


if __name__ == "__main__":
    sys.exit(main())
