#!/usr/bin/env python3
"""Exact d=3 high-charge support certificate.

This verifier uses only integer polynomial identities and Z3's exact
nonlinear-real engine.  A complex variable is represented by two real
variables; hence every UNSAT result is a statement over C, not a numerical
sample.  The predicate used below is a necessary condition for membership
in the quaternary matching-product class B_times^0([4]).

Run through verify_high_charge_q6.py together with
verify_high_charge_q6_orbits.cpp.  The C++ verifier proves
that the one-loop zero-pair equations have precisely 24 signed-linear
orbits, emits them with an explicit coordinate order, and this stage checks
their canonical-set equality with the 24 rows consumed below.  This file then
eliminates every orbit whose central support meets two different complementary
blocks.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

if not __debug__:
    raise SystemExit("exact verifier requires assertions; do not use python -O")


ZERO = -64


ORBIT_TEXT = r"""
-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,-64,-64,-64,10,11,-64,-64,-64,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,-64,-64,9,-9,11,-11,-64,-64,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,-64,-64,9,10,11,12,-64,-64,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,7,7,-64,-7,7,-64,-7,-7,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,7,7,-64,-7,11,-64,-11,-11,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,7,8,-64,-8,11,-64,-11,14,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,-64,7,8,-64,10,11,-64,13,14,-64,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,6,-6,-6,9,-64,-64,12,13,13,-13,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,6,-6,8,9,-64,-64,12,13,14,-14,-64,-64,-64,-64,-64
-64,-64,-64,-64,-64,6,7,8,9,-64,-64,12,13,14,15,-64,-64,-64,-64,-64
1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,-1
1,-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,1
1,-1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,1,1,1,1,-1
1,-1,-1,-1,-1,-1,1,8,-1,1,-1,1,13,-1,1,1,1,1,1,-1
1,-1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1
1,-1,-1,-1,-1,1,1,8,-8,1,1,12,-12,1,1,-1,-1,-1,-1,1
1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1
1,-1,-1,1,-1,1,-1,1,1,10,11,-1,-1,1,-1,1,-1,1,1,-1
1,-1,-1,1,-1,1,-1,8,-8,-1,1,12,-12,1,-1,1,-1,1,1,-1
1,-1,-1,1,-1,1,-1,8,9,-1,1,12,13,1,-1,1,-1,1,1,-1
1,-1,-1,1,-1,6,-1,1,6,-6,6,-6,-1,1,-6,1,-1,1,1,-1
1,-1,-1,1,5,-5,-1,1,-5,10,-5,10,-1,1,10,-10,-1,1,1,-1
1,-1,-1,1,5,-5,-1,1,-5,10,11,5,-1,1,5,-5,-1,1,1,-1
"""


ORBITS = [tuple(map(int, line.split(","))) for line in ORBIT_TEXT.split()]
assert len(ORBITS) == 24 and all(len(row) == 20 for row in ORBITS)


U = tuple(range(6))
TRIPLES = tuple(
    sum(1 << i for i in subset) for subset in itertools.combinations(U, 3)
)
TRIPLE_INDEX = {mask: i for i, mask in enumerate(TRIPLES)}
COMPLEMENT_INDEX = tuple(TRIPLE_INDEX[63 ^ mask] for mask in TRIPLES)
EO_WORDS = tuple(
    (x, y)
    for x in itertools.product((0, 1), repeat=2)
    for y in itertools.product((0, 1), repeat=2)
    if sum(x) + sum(y) == 2
)
PAIR_INDEX = ((0, 5), (1, 4), (2, 3))


def normalize_signed_space(row: Sequence[int]) -> tuple[int, ...]:
    """Canonical parameter names for one signed coordinate subspace."""

    components: dict[int, list[tuple[int, int]]] = {}
    for index, value in enumerate(row):
        if value != ZERO:
            components.setdefault(abs(value), []).append(
                (index, 1 if value > 0 else -1)
            )
    normalized = [ZERO] * len(row)
    for members in components.values():
        members.sort()
        anchor, anchor_sign = members[0]
        for index, sign in members:
            normalized[index] = (anchor + 1) * anchor_sign * sign
    return tuple(normalized)


def canonical_signed_space(row: Sequence[int]) -> tuple[int, ...]:
    """Canonicalize under the S_6 port action and complementation."""

    best = normalize_signed_space(row)
    for permutation in itertools.permutations(U):
        for use_complement in (False, True):
            transformed = [ZERO] * len(TRIPLES)
            for index, value in enumerate(row):
                if value == ZERO:
                    continue
                mask = 0
                for coordinate in U:
                    if (TRIPLES[index] >> coordinate) & 1:
                        mask |= 1 << permutation[coordinate]
                if use_complement:
                    mask ^= 63
                transformed[TRIPLE_INDEX[mask]] = value
            best = min(best, normalize_signed_space(transformed))
    return best


def verify_cpp_orbit_handshake(path: Path) -> None:
    """Reindex and compare the C++-generated and Python-consumed orbit sets."""

    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or lines[0] != "format=high-charge-q6-orbits-v1":
        raise AssertionError("q6 orbit handshake has an unknown format")
    if len(lines) < 3 or lines[1] != (
        "coordinate_order=increasing-integer-bitmask;bit-k-is-port-k"
    ):
        raise AssertionError("q6 orbit handshake has an unknown coordinate order")
    if not lines[2].startswith("triple_masks="):
        raise AssertionError("q6 orbit handshake is missing triple_masks")
    source_masks = tuple(map(int, lines[2].split("=", 1)[1].split(",")))
    if len(source_masks) != len(TRIPLES) or set(source_masks) != set(TRIPLES):
        raise AssertionError("q6 orbit handshake has the wrong triple coordinates")
    source_index = {mask: index for index, mask in enumerate(source_masks)}
    source_rows = [
        tuple(map(int, line.split("=", 1)[1].split(",")))
        for line in lines[3:]
        if line.startswith("orbit=")
    ]
    if len(source_rows) != 24 or any(len(row) != 20 for row in source_rows):
        raise AssertionError("q6 orbit handshake must contain 24 rows of length 20")
    reindexed = [
        tuple(row[source_index[mask]] for mask in TRIPLES) for row in source_rows
    ]
    generated = {canonical_signed_space(row) for row in reindexed}
    consumed = {canonical_signed_space(row) for row in ORBITS}
    if len(generated) != 24 or len(consumed) != 24 or generated != consumed:
        raise AssertionError(
            "q6 C++/Python canonical orbit sets disagree after coordinate reindexing"
        )
    print("PASS q6 C++/Python handshake: 24 canonical orbit sets agree")


RawCard = tuple[tuple[tuple[int, int], ...], ...]
CardMeta = tuple[tuple[int, int], tuple[int, int], tuple[int, int, int, int]]
ComplexExpr = tuple[str, str]


def generate_raw_cards() -> tuple[list[RawCard], list[CardMeta]]:
    cards: list[RawCard] = []
    metadata: list[CardMeta] = []
    seen: set[RawCard] = set()
    for first_external in itertools.combinations(U, 2):
        first_internal = tuple(i for i in U if i not in first_external)
        for second_external in itertools.combinations(U, 2):
            second_internal = tuple(i for i in U if i not in second_external)
            for second_order in itertools.permutations(second_internal):
                wiring = dict(zip(first_internal, second_order))
                entries: list[tuple[tuple[int, int], ...]] = []
                for x, y in EO_WORDS:
                    terms: list[tuple[int, int]] = []
                    for chosen in itertools.combinations(
                        first_internal, 3 - sum(x)
                    ):
                        first_mask = sum(
                            1 << first_external[k] for k in range(2) if x[k]
                        ) | sum(1 << i for i in chosen)
                        second_chosen = set(second_internal) - {
                            wiring[i] for i in chosen
                        }
                        second_mask = sum(
                            1 << second_external[k] for k in range(2) if y[k]
                        ) | sum(1 << i for i in second_chosen)
                        terms.append(
                            (
                                TRIPLE_INDEX[first_mask],
                                TRIPLE_INDEX[second_mask],
                            )
                        )
                    entries.append(tuple(terms))
                raw = tuple(entries)
                if raw not in seen:
                    seen.add(raw)
                    cards.append(raw)
                    metadata.append(
                        (first_external, second_external, second_order)
                    )
    assert len(cards) == 5400
    return cards, metadata


RAW_CARDS, RAW_CARD_METADATA = generate_raw_cards()


# Choosing the two one-positions lexicographically lists the EO words in
# reverse W4 order. The complementary-pair predicate is invariant under
# this complete reversal; two-copy arrays use W4 order directly.
LOOPS: list[tuple[tuple[int, int], ...]] = []
for i, j in itertools.combinations(U, 2):
    external = tuple(k for k in U if k not in (i, j))
    entries = []
    for pair in itertools.combinations(external, 2):
        mask = sum(1 << k for k in pair)
        entries.append(
            (TRIPLE_INDEX[mask | (1 << i)], TRIPLE_INDEX[mask | (1 << j)])
        )
    LOOPS.append(tuple(entries))


def real_add(parts: Iterable[str]) -> str:
    kept = [part for part in parts if part != "0"]
    if not kept:
        return "0"
    if len(kept) == 1:
        return kept[0]
    return "(+ " + " ".join(kept) + ")"


def real_mul(left: str, right: str) -> str:
    if left == "0" or right == "0":
        return "0"
    return f"(* {left} {right})"


def real_necessary_safe(values: Sequence[str]) -> str:
    products = [real_mul(values[a], values[b]) for a, b in PAIR_INDEX]
    clauses = [
        f"(= (= {values[a]} 0) (= {values[b]} 0))" for a, b in PAIR_INDEX
    ]
    clauses.append("(or " + " ".join(f"(= {p} 0)" for p in products) + ")")
    for i, j in itertools.combinations(range(3), 2):
        clauses.append(
            f"(or (= {products[i]} 0) (= {products[j]} 0) "
            f"(= (* {products[i]} {products[i]}) "
            f"(* {products[j]} {products[j]})))"
        )
    return "(and " + " ".join(clauses) + ")"


def orbit_roots(orbit: Sequence[int]) -> list[int]:
    return sorted({abs(value) for value in orbit if value != ZERO})


def real_q_values(orbit: Sequence[int]) -> list[str]:
    names = {root: f"x{i}" for i, root in enumerate(orbit_roots(orbit))}
    result = []
    for value in orbit:
        if value == ZERO:
            result.append("0")
        elif value > 0:
            result.append(names[abs(value)])
        else:
            result.append(f"(- {names[abs(value)]})")
    return result


def real_card_values(raw: RawCard, q: Sequence[str]) -> list[str]:
    values = []
    for entry, terms in enumerate(raw):
        summands = (["t"] if entry in (0, 5) else []) + [
            real_mul(q[i], q[j]) for i, j in terms
        ]
        values.append(real_add(summands))
    return values


def deduplicated_card_deck(orbit: Sequence[int]) -> list[RawCard]:
    q = real_q_values(orbit)
    seen: set[str] = set()
    deck: list[RawCard] = []
    for raw in RAW_CARDS:
        expression = real_necessary_safe(real_card_values(raw, q))
        if expression not in seen:
            seen.add(expression)
            deck.append(raw)
    return deck


def cneg(value: ComplexExpr) -> ComplexExpr:
    return (f"(- {value[0]})", f"(- {value[1]})")


def cadd(parts: Iterable[ComplexExpr]) -> ComplexExpr:
    parts = list(parts)
    return (
        real_add(part[0] for part in parts),
        real_add(part[1] for part in parts),
    )


def cmul(left: ComplexExpr, right: ComplexExpr) -> ComplexExpr:
    if left == ("0", "0") or right == ("0", "0"):
        return ("0", "0")
    return (
        f"(- (* {left[0]} {right[0]}) (* {left[1]} {right[1]}))",
        f"(+ (* {left[0]} {right[1]}) (* {left[1]} {right[0]}))",
    )


def czero(value: ComplexExpr) -> str:
    return f"(and (= {value[0]} 0) (= {value[1]} 0))"


def cequal(left: ComplexExpr, right: ComplexExpr) -> str:
    return f"(and (= {left[0]} {right[0]}) (= {left[1]} {right[1]}))"


def complex_necessary_safe(values: Sequence[ComplexExpr]) -> str:
    products = [cmul(values[a], values[b]) for a, b in PAIR_INDEX]
    clauses = [
        f"(= {czero(values[a])} {czero(values[b])})" for a, b in PAIR_INDEX
    ]
    clauses.append("(or " + " ".join(czero(p) for p in products) + ")")
    for i, j in itertools.combinations(range(3), 2):
        clauses.append(
            f"(or {czero(products[i])} {czero(products[j])} "
            f"{cequal(cmul(products[i], products[i]), cmul(products[j], products[j]))})"
        )
    return "(and " + " ".join(clauses) + ")"


@dataclass
class ComplexOrbit:
    q: list[ComplexExpr]
    parameters: list[ComplexExpr]


def complex_orbit(orbit: Sequence[int]) -> ComplexOrbit:
    roots = orbit_roots(orbit)
    names = {root: (f"x{i}r", f"x{i}i") for i, root in enumerate(roots)}
    q: list[ComplexExpr] = []
    for value in orbit:
        if value == ZERO:
            q.append(("0", "0"))
        elif value > 0:
            q.append(names[abs(value)])
        else:
            q.append(cneg(names[abs(value)]))
    return ComplexOrbit(q=q, parameters=[names[root] for root in roots])


def complex_card_values(
    raw: RawCard, q: Sequence[ComplexExpr]
) -> list[ComplexExpr]:
    t = ("tr", "ti")
    values = []
    for entry, terms in enumerate(raw):
        summands = ([t] if entry in (0, 5) else []) + [
            cmul(q[i], q[j]) for i, j in terms
        ]
        values.append(cadd(summands))
    return values


def complex_loop_values(
    loop: Sequence[tuple[int, int]], q: Sequence[ComplexExpr]
) -> list[ComplexExpr]:
    return [cadd((q[i], q[j])) for i, j in loop]


def target_expression(q: Sequence[ComplexExpr]) -> str:
    blocks = []
    for i, j in enumerate(COMPLEMENT_INDEX):
        if i < j:
            blocks.append(f"(or (not {czero(q[i])}) (not {czero(q[j])}))")
    return "(or " + " ".join(
        f"(and {blocks[i]} {blocks[j]})"
        for i, j in itertools.combinations(range(10), 2)
    ) + ")"


EXPECTED_Z3_VERSION = "4.14.0"
Z3: str | None = None
AUDIT_DIR: Path | None = None
AUDIT_RECORDS: list[dict] = []


def resolve_z3(requested: str | None) -> tuple[str, str]:
    candidate = requested or shutil.which("z3")
    if candidate is None:
        raise SystemExit(
            "z3 is required for the exact complex-algebra checks; "
            "install Z3 4.14.0 or pass --z3 PATH"
        )
    resolved = shutil.which(candidate) or candidate
    if not Path(resolved).is_file():
        raise SystemExit(f"z3 executable not found: {resolved}")
    process = subprocess.run(
        [resolved, "--version"],
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )
    version_text = (process.stdout or process.stderr).strip()
    match = re.search(r"\b(\d+\.\d+\.\d+)\b", version_text)
    if process.returncode != 0 or match is None:
        raise SystemExit(f"could not determine z3 version: {version_text!r}")
    if match.group(1) != EXPECTED_Z3_VERSION:
        raise SystemExit(
            f"q6 certificate pins Z3 {EXPECTED_Z3_VERSION}, "
            f"but {match.group(1)} was found at {resolved}"
        )
    return resolved, version_text


def check_unsat(lines: Sequence[str], label: str, timeout: int = 60) -> None:
    if Z3 is None:
        raise RuntimeError("internal error: z3 was not initialized")
    script = "\n".join(["(set-logic QF_NRA)", *lines, "(check-sat)"]) + "\n"
    query_name = f"query-{len(AUDIT_RECORDS) + 1:03d}.smt2"
    if AUDIT_DIR is not None:
        (AUDIT_DIR / query_name).write_text(script, encoding="utf-8", newline="\n")
    started = time.monotonic()
    process = subprocess.run(
        [Z3, "-in", f"-T:{timeout}"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout + 5,
        check=False,
    )
    status = process.stdout.splitlines()[0] if process.stdout else ""
    if AUDIT_DIR is not None:
        (AUDIT_DIR / query_name.replace('.smt2', '.stdout.txt')).write_text(
            process.stdout, encoding="utf-8", newline="\n")
        (AUDIT_DIR / query_name.replace('.smt2', '.stderr.txt')).write_text(
            process.stderr, encoding="utf-8", newline="\n")
        AUDIT_RECORDS.append({
            "file": query_name, "label": label,
            "sha256": hashlib.sha256(script.encode("utf-8")).hexdigest(),
            "result": status, "exit_code": process.returncode,
            "runtime_seconds": round(time.monotonic() - started, 6),
        })
        (AUDIT_DIR / "queries.json").write_text(
            json.dumps(AUDIT_RECORDS, indent=2) + "\n", encoding="utf-8")
    if process.returncode != 0 or status != "unsat":
        raise RuntimeError(
            f"{label}: expected unsat, got {status!r}\n{process.stdout}\n{process.stderr}"
        )


def declarations(number_of_parameters: int) -> list[str]:
    names = [f"x{i}{part}" for i in range(number_of_parameters) for part in "ri"]
    names += ["tr", "ti"]
    return [f"(declare-const {name} Real)" for name in names]


# Exact sparse-orbit cores.  Card numbers refer to the deterministic
# coefficient-expression deck constructed by deduplicated_card_deck().
SPARSE_CORES: dict[int, tuple[tuple[int, ...], tuple[int, ...]]] = {
    1: ((), ()),
    2: ((0,), (6,)),
    3: ((0,), (6,)),
    4: ((), (0,)),
    5: ((0,), (0,)),
    6: ((0,), (92,)),
    7: ((0, 1, 2), (17, 80, 92)),
    8: ((1,), (2,)),
    9: ((1, 2), (21, 208, 504)),
    10: ((1, 2), (18, 21, 208, 355, 356, 504)),
    11: ((), (0,)),
    12: ((), (0,)),
    13: ((), (0,)),
    14: ((), (600,)),
    15: ((), (0,)),
    16: (
        (0,),
        (6, 16, 17, 40, 46, 74, 76, 90, 100, 184, 190, 208,
         241, 247, 552, 1004, 1256, 3775),
    ),
    17: ((), (0,)),
    18: ((), (396,)),
    19: ((), (1, 6, 83, 101, 114, 241, 243, 246, 552, 1004,
                1252, 1253, 3775)),
}


# A few projective branches are deliberately omitted: their parameters alone
# touch only one complementary block.  The target therefore guarantees that
# at least one checked parameter is nonzero.
CHECKED_PARAMETERS: dict[int, tuple[int, ...]] = {
    7: (1, 2, 3, 4),
    9: (0, 1, 2, 4, 5),
    10: (0, 1, 3, 4, 6, 7),
}


def blocks_touched_by_parameters(
    orbit: Sequence[int], parameter_indices: set[int]
) -> int:
    roots = orbit_roots(orbit)
    selected_roots = {roots[i] for i in parameter_indices}
    touched = set()
    for i, value in enumerate(orbit):
        if value != ZERO and abs(value) in selected_roots:
            touched.add(min(i, COMPLEMENT_INDEX[i]))
    return len(touched)


def verify_sparse_orbits() -> None:
    for orbit_number in range(1, 20):
        orbit = ORBITS[orbit_number]
        current = complex_orbit(orbit)
        loop_ids, card_ids = SPARSE_CORES[orbit_number]
        deck = deduplicated_card_deck(orbit)
        if any(card_id >= len(deck) for card_id in card_ids):
            raise AssertionError(f"orbit {orbit_number}: invalid card id")

        base = declarations(len(current.parameters))
        base += [
            f"(assert (not {czero(('tr', 'ti'))}))",
            f"(assert {target_expression(current.q)})",
        ]
        for loop_id in loop_ids:
            base.append(
                f"(assert {complex_necessary_safe(complex_loop_values(LOOPS[loop_id], current.q))})"
            )
        for card_id in card_ids:
            base.append(
                f"(assert {complex_necessary_safe(complex_card_values(deck[card_id], current.q))})"
            )

        checked = CHECKED_PARAMETERS.get(
            orbit_number, tuple(range(len(current.parameters)))
        )
        skipped = set(range(len(current.parameters))) - set(checked)
        if skipped and blocks_touched_by_parameters(orbit, skipped) > 1:
            raise AssertionError(
                f"orbit {orbit_number}: skipped parameters can satisfy the target"
            )

        for parameter in checked:
            real, imaginary = current.parameters[parameter]
            check_unsat(
                [
                    *base,
                    f"(assert (= {real} 1))",
                    f"(assert (= {imaginary} 0))",
                ],
                f"orbit {orbit_number}, parameter {parameter}",
            )
        print(f"PASS sparse orbit {orbit_number}")


def raw_card_from_meta(meta: CardMeta) -> RawCard:
    try:
        return RAW_CARDS[RAW_CARD_METADATA.index(meta)]
    except ValueError as error:
        raise AssertionError(f"unknown card metadata {meta}") from error


def implication_check(
    orbit_number: int,
    card_meta: Sequence[CardMeta],
    assumptions: Sequence[str],
    label: str,
) -> None:
    current = complex_orbit(ORBITS[orbit_number])
    lines = declarations(len(current.parameters))
    lines.append(f"(assert (not {czero(('tr', 'ti'))}))")
    lines.extend(assumptions)
    for meta in card_meta:
        raw = raw_card_from_meta(meta)
        lines.append(
            f"(assert {complex_necessary_safe(complex_card_values(raw, current.q))})"
        )
    check_unsat(lines, label)


def nonzero(parameter: int) -> str:
    return f"(assert (not {czero((f'x{parameter}r', f'x{parameter}i'))}))"


def zero(parameter: int) -> str:
    return f"(assert {czero((f'x{parameter}r', f'x{parameter}i'))})"


def verify_dense_orbits() -> None:
    identity_03: CardMeta = ((0, 3), (0, 3), (1, 2, 4, 5))

    # Orbit 20: this card is (E,2a^2,-6a^2,-6a^2,2a^2,E),
    # so it forces the 16-entry component a=x0 to vanish.  The residual
    # four singleton coordinates are a port permutation of sparse orbit 3.
    implication_check(20, [identity_03], [nonzero(0)], "orbit 20, dense component")
    residual_positions = {
        i for i, value in enumerate(ORBITS[20]) if value not in (ZERO, 1, -1)
    }
    orbit3_positions = {i for i, value in enumerate(ORBITS[3]) if value != ZERO}
    permutation = (0, 1, 2, 5, 3, 4)
    transformed = set()
    for index in residual_positions:
        mask = 0
        for coordinate in U:
            if (TRIPLES[index] >> coordinate) & 1:
                mask |= 1 << permutation[coordinate]
        transformed.add(TRIPLE_INDEX[mask])
    if transformed != orbit3_positions:
        raise AssertionError("orbit 20 residual is not orbit 3")
    print("PASS dense orbit 20 -> sparse orbit 3")

    # Orbit 21: the two cards are respectively
    # (t,2b^2,-2b^2,-2b^2,2b^2,t) and the same expression with a.
    card21_a: CardMeta = ((0, 1), (0, 1), (3, 4, 5, 2))
    card21_b: CardMeta = ((0, 4), (0, 4), (2, 3, 5, 1))
    implication_check(21, [card21_a], [nonzero(1)], "orbit 21, parameter 1")
    implication_check(21, [card21_b], [nonzero(0)], "orbit 21, parameter 0")
    print("PASS dense orbit 21")

    # Orbit 22: the first two cards force x1*x2=0 and x1^2+x2^2=0,
    # hence x1=x2=0 over C.  The identity_03 card then forces x0=0.
    card22_product: CardMeta = ((0, 1), (0, 1), (2, 3, 5, 4))
    card22_sum: CardMeta = ((0, 1), (0, 1), (4, 5, 3, 2))
    implication_check(
        22,
        [card22_product, card22_sum],
        [
            f"(assert (or (not {czero(('x1r', 'x1i'))}) "
            f"(not {czero(('x2r', 'x2i'))})))"
        ],
        "orbit 22, secondary components",
    )
    implication_check(
        22,
        [identity_03],
        [zero(1), zero(2), nonzero(0)],
        "orbit 22, principal component",
    )
    print("PASS dense orbit 22")

    # Orbit 23 has a three-card hand certificate.  The first card forces
    # b=x1=0.  With b=0, the second gives a!=0 => c=d=0.  The third then
    # rules out a!=0.  The surviving c,d entries are complementary triples.
    card23_first: CardMeta = ((0, 1), (0, 1), (3, 4, 5, 2))
    card23_second: CardMeta = ((0, 1), (1, 4), (2, 3, 0, 5))
    card23_third: CardMeta = ((0, 3), (0, 3), (5, 1, 2, 4))
    current23 = complex_orbit(ORBITS[23])
    common23 = declarations(4) + [f"(assert (not {czero(('tr', 'ti'))}))"]
    first_values = complex_card_values(
        raw_card_from_meta(card23_first), current23.q
    )

    # If b!=0, the first middle pair (entries 1 and 4) is active.
    # Since the endpoint pair is active as well, safety would force entries
    # 2 and 3 both to vanish.  Their exact formulas make that impossible.
    check_unsat(
        [
            *common23,
            nonzero(1),
            f"(assert (or {czero(first_values[1])} {czero(first_values[4])}))",
        ],
        "orbit 23, first-card middle pair is active",
    )
    check_unsat(
        [
            *common23,
            nonzero(1),
            f"(assert {czero(first_values[2])})",
            f"(assert {czero(first_values[3])})",
        ],
        "orbit 23, first card forces b=0",
    )

    second_values = complex_card_values(
        raw_card_from_meta(card23_second), current23.q
    )
    cd = cmul(current23.parameters[2], current23.parameters[3])
    c_or_d_nonzero = (
        f"(assert (or (not {czero(current23.parameters[2])}) "
        f"(not {czero(current23.parameters[3])})))"
    )

    # With b=0 the second card is (t,cd,ad,-ac,cd,t).  If cd!=0,
    # endpoint and first-middle pairs are active, so the second-middle pair
    # would have to vanish.  If cd=0, complementary-pair closure of that
    # second-middle pair is already impossible unless c=d=0.
    check_unsat(
        [
            *common23,
            zero(1),
            nonzero(0),
            f"(assert (not {czero(cd)}))",
            f"(assert {czero(second_values[2])})",
            f"(assert {czero(second_values[3])})",
        ],
        "orbit 23, cd!=0 branch",
    )
    check_unsat(
        [
            *common23,
            zero(1),
            nonzero(0),
            c_or_d_nonzero,
            f"(assert {czero(cd)})",
            f"(assert (= {czero(second_values[2])} {czero(second_values[3])}))",
        ],
        "orbit 23, cd=0 branch",
    )

    third_values = complex_card_values(
        raw_card_from_meta(card23_third), current23.q
    )
    third_products = [
        cmul(third_values[a], third_values[b]) for a, b in PAIR_INDEX
    ]
    # Under b=c=d=0 this card is
    # (t,3a^2,-3a^2,-3a^2,3a^2,t).  Safety requires at least one
    # complementary pair product to vanish, contradicting a*t!=0.
    check_unsat(
        [
            *common23,
            zero(1),
            zero(2),
            zero(3),
            nonzero(0),
            "(assert (or "
            + " ".join(czero(product) for product in third_products)
            + "))",
        ],
        "orbit 23, principal component",
    )
    roots = orbit_roots(ORBITS[23])
    residual_root_set = {roots[2], roots[3]}
    residual_indices = [
        i
        for i, value in enumerate(ORBITS[23])
        if value != ZERO and abs(value) in residual_root_set
    ]
    if len(residual_indices) != 2 or COMPLEMENT_INDEX[residual_indices[0]] != residual_indices[1]:
        raise AssertionError("orbit 23 residual is not one complementary pair")
    print("PASS dense orbit 23 -> one complementary pair")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replay the exact q6 algebra stage after the C++ orbit audit."
    )
    parser.add_argument(
        "--cpp-orbits",
        required=True,
        type=Path,
        help="orbit handshake emitted by verify_high_charge_q6_orbits.cpp",
    )
    parser.add_argument(
        "--z3",
        help="path to the pinned Z3 4.14.0 executable (otherwise use PATH)",
    )
    parser.add_argument("--audit-dir", type=Path,
                        help="retain every exact SMT input and solver response")
    args = parser.parse_args()

    verify_cpp_orbit_handshake(args.cpp_orbits)
    global Z3, AUDIT_DIR
    Z3, version_text = resolve_z3(args.z3)
    AUDIT_DIR = args.audit_dir
    if AUDIT_DIR is not None:
        AUDIT_DIR.mkdir(parents=True, exist_ok=True)
        if (AUDIT_DIR / "queries.json").exists():
            raise RuntimeError("audit directory already contains a run; use a new directory")
        metadata = {
            "format": "paper2-q6-exact-inputs-v1",
            "z3_version": version_text,
            "z3_executable_sha256": hashlib.sha256(Path(Z3).read_bytes()).hexdigest(),
            "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "port_labels": "0,...,5; bit k is port k",
            "triple_masks_in_coefficient_order": TRIPLES,
            "eo_words": EO_WORDS,
            "orbit_rows": ORBITS,
            "zero_sentinel": ZERO,
            "orbit_parameter_roots": [orbit_roots(row) for row in ORBITS],
            "sparse_cores": SPARSE_CORES,
            "checked_parameter_overrides": CHECKED_PARAMETERS,
            "dense_branch_rules": "verify_dense_orbits in the hashed generator",
            "proof_objects": False,
            "trust_boundary": "Z3 nonlinear-real UNSAT answers are trusted, not independently proof-checked",
        }
        (AUDIT_DIR / "encoding.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Z3: {version_text} ({Z3})")
    verify_sparse_orbits()
    verify_dense_orbits()
    if AUDIT_DIR is not None:
        print(f"archived_exact_queries: {len(AUDIT_RECORDS)} ({AUDIT_DIR})")
    print("PASS: every d=3 joint-safe branch is contained in one complementary pair")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, AssertionError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
