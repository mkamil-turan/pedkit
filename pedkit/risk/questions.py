"""Answer the ``##question`` lines of a PEDX file.

    ##question=<ID=Q1,type=mode,disease=D1>                          -> mode inference + diagnosis
    ##question=<ID=Q2,type=risk,individual=II-2,disease=D1,event=carrier> -> Bayes risk
    (for a disease with mode=MF a risk question is answered by the
    liability-threshold model, ``pedkit.mf``)
"""

from __future__ import annotations

from dataclasses import dataclass

from ..inference.exact import infer_mode_exact
from ..inference.montecarlo import infer_mode
from ..model import Mode, Pedigree, Question
from .bayes import bayes_risk


@dataclass
class QuestionAnswer:
    question: Question
    result: object | None   # ModeInference, RiskResult or MFRiskResult
    error: str | None = None

    def report(self) -> str:
        q = self.question
        title = f"Question {q.id or '?'} (type: {q.type})"
        if self.error:
            return f"{title}\n  Could not be answered: {self.error}"
        return f"{title}\n\n{self.result.report()}"


def answer_question(ped: Pedigree, question: Question, *, n: int = 1_000_000,
                    seed: int | None = None, method: str = "auto",
                    proposal_freq: float | None = None) -> QuestionAnswer:
    try:
        if question.type == "mode":
            if method == "monte-carlo":
                result = infer_mode(ped, question.disease, n=n, seed=seed,
                                    proposal_freq=proposal_freq)
            else:
                result = infer_mode_exact(ped, question.disease)
            return QuestionAnswer(question, result)
        if question.type == "risk":
            if not question.individual:
                return QuestionAnswer(question, None, "the risk question has no 'individual' field")
            d = ped.disease(question.disease) if question.disease else ped.primary_disease
            if d is not None and d.mode == Mode.MF:
                from ..mf import mf_risk
                return QuestionAnswer(question, mf_risk(ped, question.individual, d, seed=seed))
            result = bayes_risk(ped, question.individual, question.event or "carrier",
                                question.disease, method=method, n=n, seed=seed,
                                proposal_freq=proposal_freq)
            return QuestionAnswer(question, result)
        return QuestionAnswer(question, None, f"unknown question type {question.type!r} "
                                              "(expected 'mode' or 'risk')")
    except (KeyError, ValueError) as e:
        return QuestionAnswer(question, None, str(e))


def answer_questions(ped: Pedigree, **options) -> list[QuestionAnswer]:
    return [answer_question(ped, q, **options) for q in ped.questions]
