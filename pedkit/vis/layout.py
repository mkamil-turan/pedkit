"""Pedigree layout: assigns each individual a generation and a horizontal position.

A simple algorithm for pedigrees of 3-5 generations and fewer than about 30
individuals:

1. Generation: one more than the deepest parent; spouses are pulled to the
   same generation; spouses whose parents are not drawn join the generation
   of their partner.
2. Sibships (same mother and father) are ordered by ``order`` or file order,
   and twins are kept next to each other.
3. A spouse marrying into the family is placed to the right of the
   individual from the family.
4. Generations are placed top-down, left to right; a sibship is centred
   under its parents and shifted right on overlap. A bottom-up pass then
   tries to centre parents above their children.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..model import Individual, Pedigree

H_SPACING = 80.0      # distance between neighbouring symbol centres
BLOCK_GAP = 20.0      # extra space between different sibship blocks
FAMILY_GAP = 120.0    # space between families when the file holds several

Key = tuple[str, str]  # (family, individual)


def key(ind: Individual) -> Key:
    return (ind.family, ind.id)


@dataclass
class Sibship:
    parents: list[Individual]
    children: list[Individual]


@dataclass
class Layout:
    x: dict[Key, float] = field(default_factory=dict)
    generation: dict[Key, int] = field(default_factory=dict)
    individuals: list[Individual] = field(default_factory=list)  # to draw (no duplicates)
    couples: list[tuple[Individual, Individual]] = field(default_factory=list)
    sibships: list[Sibship] = field(default_factory=list)

    @property
    def n_generations(self) -> int:
        return max(self.generation.values(), default=-1) + 1


def compute_layout(ped: Pedigree) -> Layout:
    result = Layout()
    left = 0.0
    for family in ped.families:
        part = _FamilyLayout(ped, family).run()
        if not part.x:
            continue
        shift = left - min(part.x.values())
        for k in part.x:
            part.x[k] += shift
        result.x.update(part.x)
        result.generation.update(part.generation)
        result.individuals += part.individuals
        result.couples += part.couples
        result.sibships += part.sibships
        left = max(part.x.values()) + FAMILY_GAP
    return result


@dataclass
class _Block:
    items: list[str]                   # individual IDs, left to right
    siblings: list[str]                # sibship members among items
    parents: tuple[str, ...] = ()      # parents of the sibship


class _FamilyLayout:
    def __init__(self, ped: Pedigree, family: str):
        self.ped = ped
        self.family = family
        self.members: dict[str, Individual] = {}
        for ind in ped.individuals:
            if ind.family == family and ind.id not in self.members:  # duplicates skipped
                self.members[ind.id] = ind
        self.file_order = {iid: i for i, iid in enumerate(self.members)}

    # -- relationships ----------------------------------------------------------

    def parents(self, iid: str) -> tuple[str, ...]:
        ind = self.members[iid]
        return tuple(p for p in (ind.father, ind.mother) if p in self.members and p != iid)

    def is_founder(self, iid: str) -> bool:
        return not self.parents(iid)

    def _build_relations(self) -> None:
        self.spouses: dict[str, list[str]] = {iid: [] for iid in self.members}
        self.couples: list[tuple[str, str]] = []

        def add_couple(a: str, b: str) -> None:
            if a == b or b in self.spouses[a]:
                return
            self.spouses[a].append(b)
            self.spouses[b].append(a)
            self.couples.append((a, b))

        self.sibship: dict[tuple[str, ...], list[str]] = {}
        for iid in self.members:
            ps = self.parents(iid)
            if ps:
                self.sibship.setdefault(ps, []).append(iid)
            if len(ps) == 2:
                add_couple(*ps)
        for u in self.ped.unions:
            if u.p1 in self.members and u.p2 in self.members:
                add_couple(u.p1, u.p2)
        for ps, kids in self.sibship.items():
            self.sibship[ps] = self._order_siblings(kids)

    def _order_siblings(self, kids: list[str]) -> list[str]:
        inf = float("inf")
        ordered = sorted(kids, key=lambda c: (
            self.members[c].order if self.members[c].order is not None else inf,
            self.file_order[c]))
        # gather members of a twin group next to the group's first member
        result: list[str] = []
        for c in ordered:
            if c in result:
                continue
            result.append(c)
            twin = self.members[c].twin
            if twin:
                result += [d for d in ordered if d not in result and self.members[d].twin
                           and self.members[d].twin.group == twin.group]
        return result

    # -- 1. generation ---------------------------------------------------------

    def _assign_generations(self) -> dict[str, int]:
        g = {iid: 0 for iid in self.members}
        for _ in range(2 * len(self.members) + 2):  # bounded even if the pedigree loops
            changed = False

            def raise_to(iid: str, value: int) -> None:
                nonlocal changed
                if g[iid] < value:
                    g[iid] = value
                    changed = True

            for iid in self.members:
                ps = self.parents(iid)
                if ps:
                    raise_to(iid, max(g[p] for p in ps) + 1)
            for a, b in self.couples:
                raise_to(a, g[b])
                raise_to(b, g[a])
            # parents who are both founders move down to just above their children
            for ps, kids in self.sibship.items():
                if all(self.is_founder(p) for p in ps):
                    for p in ps:
                        raise_to(p, min(g[c] for c in kids) - 1)
            if not changed:
                break
        return g

    # -- 2-3. ordering ---------------------------------------------------------

    def _order(self, g: dict[str, int]) -> dict[int, list[_Block]]:
        placed: set[str] = set()
        rows: dict[int, list[_Block]] = {k: [] for k in sorted(set(g.values()))}

        def unit(iid: str) -> list[str]:
            """The individual plus the married-in spouse(s) to the right."""
            items = [iid]
            placed.add(iid)
            for s in self.spouses[iid]:
                if s not in placed and self.is_founder(s) and g[s] == g[iid]:
                    items.append(s)
                    placed.add(s)
            return items

        def is_root(iid: str) -> bool:
            # a founder, unless it has a spouse from the family (then it goes next to them)
            return self.is_founder(iid) and not any(not self.is_founder(s)
                                                    for s in self.spouses[iid])

        done: set[tuple[str, ...]] = set()
        for gen in sorted(rows):
            for iid in self.members:
                if g[iid] == gen and iid not in placed and is_root(iid):
                    rows[gen].append(_Block(unit(iid), [iid]))
            for block in list(rows[gen]):
                for p in block.items:
                    for ps, kids in self.sibship.items():
                        if p not in ps or ps in done:
                            continue
                        done.add(ps)
                        pending = [c for c in kids if c not in placed]
                        for cg in sorted({g[c] for c in pending}):
                            group = [c for c in pending if g[c] == cg]
                            items: list[str] = []
                            for c in group:
                                if c not in placed:
                                    items += unit(c)
                            rows[cg].append(_Block(items, group, ps))
        for iid in self.members:  # anything left unplaced (e.g. loops)
            if iid not in placed:
                rows[g[iid]].append(_Block(unit(iid), [iid]))
        return rows

    # -- 4. horizontal position -----------------------------------------------

    def _horizontal(self, rows: dict[int, list[_Block]]) -> dict[str, float]:
        x: dict[str, float] = {}
        for gen in sorted(rows):
            cursor: float | None = None
            for block in rows[gen]:
                wanted = cursor if cursor is not None else 0.0
                placed_parents = [p for p in block.parents if p in x]
                if placed_parents:
                    mid = sum(x[p] for p in placed_parents) / len(placed_parents)
                    i0 = block.items.index(block.siblings[0])
                    i1 = block.items.index(block.siblings[-1])
                    wanted = mid - (i0 + i1) / 2 * H_SPACING
                left = wanted if cursor is None else max(cursor, wanted)
                for i, iid in enumerate(block.items):
                    x[iid] = left + i * H_SPACING
                cursor = x[block.items[-1]] + H_SPACING + BLOCK_GAP

        # bottom-up: centre parents above their children
        for gen in sorted(rows, reverse=True):
            blocks = rows[gen]
            for i in reversed(range(len(blocks))):
                block = blocks[i]
                shifts = []
                for ps, kids in self.sibship.items():
                    if any(p in block.items for p in ps) and all(c in x for c in kids):
                        kids_mid = (x[kids[0]] + x[kids[-1]]) / 2
                        shifts.append(kids_mid - sum(x[p] for p in ps) / len(ps))
                if not shifts:
                    continue
                d = sum(shifts) / len(shifts)
                if i + 1 < len(blocks):
                    d = min(d, x[blocks[i + 1].items[0]] - H_SPACING - BLOCK_GAP
                            - x[block.items[-1]])
                if i > 0:
                    d = max(d, x[blocks[i - 1].items[-1]] + H_SPACING + BLOCK_GAP
                            - x[block.items[0]])
                for iid in block.items:
                    x[iid] += d
        return x

    # -- main -----------------------------------------------------------------

    def run(self) -> Layout:
        self._build_relations()
        g = self._assign_generations()
        if g:
            lowest = min(g.values())
            g = {k: v - lowest for k, v in g.items()}
        rows = self._order(g)
        x = self._horizontal(rows)

        m = self.members
        return Layout(
            x={(self.family, k): v for k, v in x.items()},
            generation={(self.family, k): v for k, v in g.items()},
            individuals=list(m.values()),
            couples=[(m[a], m[b]) for a, b in self.couples],
            sibships=[Sibship([m[p] for p in ps], [m[c] for c in cs])
                      for ps, cs in self.sibship.items()],
        )
