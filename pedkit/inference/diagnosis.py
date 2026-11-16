"""Hierarchical inheritance diagnosis from the mode posteriors.

    Inheritance
    ├─ Autosomal ──────── dominant · recessive · codominant
    ├─ Gonosomal
    │   ├─ X-linked ───── dominant · recessive
    │   └─ Y-linked
    └─ Mitochondrial (maternal)

The probability of a node is the sum of the posteriors of the modes below
it; the conditional probability is the node's share of its parent. At every
level the diagnosis is the most probable branch, without any threshold
(equally probable branches are reported together). Only the modes tested
in the inference appear in the tree.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..model import Mode
from .result import ModeInference

MODE_NAME = {
    Mode.AD: "autosomal dominant", Mode.AR: "autosomal recessive",
    Mode.ACD: "autosomal codominant", Mode.XD: "X-linked dominant",
    Mode.XR: "X-linked recessive", Mode.Y: "Y-linked",
    Mode.MT: "mitochondrial (maternal)", Mode.MF: "multifactorial",
}

# (name, children) or a Mode as a leaf
TREE = ("Inheritance", [
    ("Autosomal", [Mode.AD, Mode.AR, Mode.ACD]),
    ("Gonosomal", [
        ("X-linked", [Mode.XD, Mode.XR]),
        ("Y-linked", [Mode.Y]),
    ]),
    ("Mitochondrial", [Mode.MT]),
])


@dataclass
class DiagnosisNode:
    name: str
    mode: Mode | None                     # set for a leaf
    probability: float | None             # P(node | data)
    conditional: float | None             # P(node | parent node, data)
    prior: float
    children: list[DiagnosisNode] = field(default_factory=list)

    @property
    def most_likely_children(self) -> list[DiagnosisNode]:
        """Child node(s) with the highest probability; empty if all are zero."""
        values = [c.probability or 0.0 for c in self.children]
        best = max(values, default=0.0)
        if best <= 0:
            return []
        return [c for c, v in zip(self.children, values) if abs(v - best) <= 1e-12]


@dataclass
class HierarchicalDiagnosis:
    root: DiagnosisNode
    path: list[list[DiagnosisNode]]  # most likely node(s) at each level

    @property
    def diagnosis(self) -> str:
        if not self.path:
            return "No mode is compatible with the observations; no diagnosis."
        last = self.path[-1]
        if len(last) == 1 and last[0].mode is not None:
            return MODE_NAME[last[0].mode]
        return " / ".join(d.name if d.mode is None else MODE_NAME[d.mode] for d in last) + \
            " (equally probable)"

    def report(self) -> str:
        lines = ["Hierarchical diagnosis", ""]
        _write_tree(self.root, "", True, lines, root=True)
        lines.append("")
        for i, level in enumerate(self.path, start=1):
            names = " / ".join(MODE_NAME[d.mode] if d.mode else d.name for d in level)
            lines.append(f"  Level {i}: {names}  (P = {_fmt(level[0].probability)}, within "
                         f"parent {_fmt(level[0].conditional)})")
        lines += ["", f"Diagnosis: {self.diagnosis}"]
        return "\n".join(lines)


def hierarchical_diagnosis(c: ModeInference) -> HierarchicalDiagnosis:
    tested = set(c.results)
    all_zero = c.most_likely is None

    def build(node, parent_p: float | None) -> DiagnosisNode | None:
        if isinstance(node, Mode):
            if node not in tested:
                return None
            r = c.results[node]
            p = None if all_zero else (r.posterior or 0.0)
            return DiagnosisNode(MODE_NAME[node], node, p, _ratio(p, parent_p), r.prior)
        name, kids = node
        leaves = [m for m in _leaves(node) if m in tested]
        if not leaves:
            return None
        p = None if all_zero else sum(c.results[m].posterior or 0.0 for m in leaves)
        prior = sum(c.results[m].prior for m in leaves)
        d = DiagnosisNode(name, None, p, _ratio(p, parent_p), prior)
        d.children = [k for k in (build(x, p) for x in kids) if k is not None]
        # collapse an intermediate node with a single leaf (e.g. "Mitochondrial" > MT)
        if len(d.children) == 1 and d.children[0].mode is not None and not d.children[0].children:
            d.mode = d.children[0].mode
            d.name = MODE_NAME[d.mode][0].upper() + MODE_NAME[d.mode][1:]
            d.children = []
        return d

    root = build(TREE, None)
    root.conditional = None
    path: list[list[DiagnosisNode]] = []
    node = root
    while node.children:
        best = node.most_likely_children
        if not best:
            break
        path.append(best)
        if len(best) > 1:
            break
        node = best[0]
    return HierarchicalDiagnosis(root, path)


def _leaves(node) -> list[Mode]:
    if isinstance(node, Mode):
        return [node]
    return [m for x in node[1] for m in _leaves(x)]


def _ratio(p: float | None, parent: float | None) -> float | None:
    if p is None:
        return None
    if parent is None:
        return p
    return p / parent if parent > 0 else 0.0


def _fmt(x: float | None) -> str:
    if x is None:
        return "—"
    if x == 0:
        return "0"
    if x < 1e-3:
        return f"{x:.2e}"
    return f"{x:.4f}".rstrip("0").rstrip(".")


def _write_tree(d: DiagnosisNode, prefix: str, last: bool, lines: list[str], root: bool = False):
    info = f"P = {_fmt(d.probability)}"
    if not root:
        info += f" · within parent {_fmt(d.conditional)} · prior {_fmt(d.prior)}"
    if root:
        lines.append(f"  {d.name}  ({info})")
        new_prefix = "  "
    else:
        lines.append(f"  {prefix}{'└─ ' if last else '├─ '}{d.name}  ({info})")
        new_prefix = prefix + ("   " if last else "│  ")
    for i, k in enumerate(d.children):
        _write_tree(k, new_prefix, i == len(d.children) - 1, lines)
