#!/usr/bin/env python3
"""Exact replay for the odd-dihedral arity-six finite claims.

The script deliberately starts from the defining q4 normal forms and from
the definition of a phase card.  It does not read a payload, a list of
survivors, or a checksum.  Arithmetic is in Q(w), w^2+w+1=0.

The support audit is independent of coefficients.  The coefficient audit
constructs the q4 atlas and all of its projective three-secants exactly.
"""

from __future__ import annotations

import argparse
import itertools
import time
from collections import Counter, defaultdict
from fractions import Fraction

if not __debug__:
    raise SystemExit("exact verifier requires assertions; do not use python -O")


class Qw:
    """Element a+b*w of Q(w), where w^2=-w-1."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, other):
        other = as_qw(other)
        return Qw(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Qw(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-as_qw(other))

    def __rsub__(self, other):
        return as_qw(other) - self

    def __mul__(self, other):
        other = as_qw(other)
        # (a+bw)(c+dw)=(ac-bd)+(ad+bc-bd)w.
        return Qw(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    __rmul__ = __mul__

    def inverse(self):
        if not self:
            raise ZeroDivisionError
        # Norm(a+bw)=(a+bw)(a+bw^2)=a^2-ab+b^2.
        norm = self.a * self.a - self.a * self.b + self.b * self.b
        return Qw((self.a - self.b) / norm, -self.b / norm)

    def __truediv__(self, other):
        return self * as_qw(other).inverse()

    def __rtruediv__(self, other):
        return as_qw(other) / self

    def __pow__(self, exponent):
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        ans, base = ONE, self
        while exponent:
            if exponent & 1:
                ans = ans * base
            base = base * base
            exponent >>= 1
        return ans

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, other):
        try:
            other = as_qw(other)
        except TypeError:
            return False
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f"Qw({self.a!r},{self.b!r})"

    def compact(self):
        return (self.a.numerator, self.a.denominator,
                self.b.numerator, self.b.denominator)


def as_qw(value):
    return value if isinstance(value, Qw) else Qw(value)


ZERO = Qw(0)
ONE = Qw(1)
W = Qw(0, 1)
W2 = W * W
MU3 = (ONE, W, W2)


def vadd(u, v):
    return tuple(a + b for a, b in zip(u, v))


def vsub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def vscale(a, v):
    return tuple(a * x for x in v)


def first_nonzero(v):
    return next((i for i, x in enumerate(v) if x), None)


def projective(v):
    pivot = first_nonzero(v)
    if pivot is None:
        raise ValueError("zero vector has no projective normalization")
    return vscale(ONE / v[pivot], v)


def line_key(v, w):
    """Canonical RREF of the 2-space spanned by independent v,w."""
    rows = [list(v), list(w)]
    pivot1 = next(i for i in range(len(v)) if rows[0][i] or rows[1][i])
    if not rows[0][pivot1]:
        rows.reverse()
    z = ONE / rows[0][pivot1]
    rows[0] = [z * x for x in rows[0]]
    c = rows[1][pivot1]
    rows[1] = [x - c * y for x, y in zip(rows[1], rows[0])]
    pivot2 = first_nonzero(rows[1])
    if pivot2 is None:
        raise ValueError("dependent vectors")
    z = ONE / rows[1][pivot2]
    rows[1] = [z * x for x in rows[1]]
    c = rows[0][pivot2]
    rows[0] = [x - c * y for x, y in zip(rows[0], rows[1])]
    return tuple(rows[0]), tuple(rows[1])


def bit_tuples(n):
    return tuple(itertools.product((0, 1), repeat=n))


BITS4 = bit_tuples(4)
BITS6 = bit_tuples(6)
INDEX4 = {x: i for i, x in enumerate(BITS4)}
INDEX6 = {x: i for i, x in enumerate(BITS6)}


def complement(x):
    return tuple(1 - b for b in x)


def complement_blocks(n, parity):
    words = [x for x in bit_tuples(n) if sum(x) % 2 == parity]
    blocks = []
    seen = set()
    for x in words:
        if x in seen:
            continue
        y = complement(x)
        seen.add(x)
        seen.add(y)
        blocks.append((x, y))  # lexicographic orientation
    return tuple(blocks)


BLOCKS4 = {parity: complement_blocks(4, parity) for parity in (0, 1)}


def q4_atlas(parity):
    """Definition-level projective atlas from the one/two-block forms."""
    points = {}
    blocks = BLOCKS4[parity]
    for p, (x, xb) in enumerate(blocks):
        for r, delta in enumerate(MU3):
            v = [ZERO] * 16
            v[INDEX4[x]], v[INDEX4[xb]] = ONE, delta
            points[projective(v)] = (1, (p,), (r,))
    for p, q in itertools.combinations(range(4), 2):
        x, xb = blocks[p]
        y, yb = blocks[q]
        for r, delta in enumerate(MU3):
            for t, gamma in enumerate(MU3):
                v = [ZERO] * 16
                v[INDEX4[x]] = ONE
                v[INDEX4[xb]] = delta
                v[INDEX4[y]] = -gamma
                v[INDEX4[yb]] = -delta / gamma
                points[projective(v)] = (2, (p, q), (r, t))
    assert len(points) == 66
    return points


def q4_line_audit():
    results = {}
    for parity in (0, 1):
        atlas = q4_atlas(parity)
        points = tuple(atlas)
        lines = defaultdict(set)
        for i, j in itertools.combinations(range(len(points)), 2):
            lines[line_key(points[i], points[j])].update((i, j))
        rich = {key: ids for key, ids in lines.items() if len(ids) >= 3}
        assert all(len(ids) == 3 for ids in rich.values())
        histogram = Counter(
            tuple(sorted(atlas[points[i]][0] for i in ids))
            for ids in rich.values()
        )
        # Split 2/2/2 into fixed two-block ruling versus three-block triangle.
        geometry = Counter()
        for ids in rich.values():
            types = tuple(sorted(atlas[points[i]][0] for i in ids))
            block_union = set()
            block_sets = []
            for i in ids:
                bs = set(atlas[points[i]][1])
                block_sets.append(bs)
                block_union |= bs
            if types == (1, 1, 1):
                geometry["one-block"] += 1
            elif types == (1, 1, 2):
                geometry["secant"] += 1
            elif types == (2, 2, 2) and len(block_union) == 2:
                geometry["ruling"] += 1
            elif types == (2, 2, 2) and len(block_union) == 3:
                geometry["triangle"] += 1
            else:
                raise AssertionError((types, block_sets))
        assert geometry == Counter(
            {"one-block": 4, "ruling": 36, "secant": 54,
             "triangle": 108}
        )
        assert len(rich) == 202
        results[parity] = (atlas, rich, geometry)
    return results


def dependence(v0, v1, v2):
    """Return nonzero r_i with r0*v0+r1*v1+r2*v2=0."""
    for i, j in itertools.combinations(range(len(v0)), 2):
        determinant = v0[i] * v1[j] - v0[j] * v1[i]
        if determinant:
            # v2 = alpha*v0 + beta*v1.
            alpha = (v2[i] * v1[j] - v2[j] * v1[i]) / determinant
            beta = (v0[i] * v2[j] - v0[j] * v2[i]) / determinant
            assert v2 == vadd(vscale(alpha, v0), vscale(beta, v1))
            relation = (-alpha, -beta, ONE)
            assert all(relation)
            return relation
    raise AssertionError("first two projective points were dependent")


def labeled_pencils(atlas, rich, geometry_filter=None):
    """Generate (A,B,geometry) modulo one common nonzero scalar.

    The values at 1,w,w^2 are not chosen from a payload.  Their relative
    scalars are the unique ones imposed by the affine-pencil circuit.
    """
    points = tuple(atlas)
    weights = (W - W2, W2 - ONE, ONE - W)
    answer = {}
    for key, ids in rich.items():
        ids = tuple(ids)
        kinds = tuple(sorted(atlas[points[i]][0] for i in ids))
        union = set().union(*(set(atlas[points[i]][1]) for i in ids))
        if kinds == (1, 1, 1):
            geometry = "one-block"
        elif kinds == (1, 1, 2):
            geometry = "secant"
        elif len(union) == 2:
            geometry = "ruling"
        else:
            geometry = "triangle"
        if geometry_filter is not None and geometry not in geometry_filter:
            continue
        for ordered in itertools.permutations(ids):
            vectors = tuple(points[i] for i in ordered)
            relation = dependence(*vectors)
            q = tuple(
                vscale(relation[i] / weights[i], vectors[i])
                for i in range(3)
            )
            # q_1=A-B and q_w=A-wB.
            B = vscale(ONE / (W - ONE), vsub(q[0], q[1]))
            A = vadd(q[0], B)
            assert q[2] == vsub(A, vscale(W2, B))
            joint_pivot = first_nonzero(A + B)
            scale = ONE / (A + B)[joint_pivot]
            A, B = vscale(scale, A), vscale(scale, B)
            answer[(A, B)] = geometry
    return answer


def atlas_contains(v, atlas_set):
    return first_nonzero(v) is None or projective(v) in atlas_set


def block_vector(parity, values):
    """Create a q4 vector from four ordered complement-block value pairs."""
    assert len(values) == 4
    v = [ZERO] * 16
    for block, pair_values in zip(BLOCKS4[parity], values):
        for x, value in zip(block, pair_values):
            v[INDEX4[x]] = value
    return tuple(v)


def switching_representatives():
    secant_A = block_vector(0, ((ONE, ONE), (W, W), (ZERO, ZERO), (ZERO, ZERO)))
    secant_B = block_vector(0, ((W2, W2), (W, W), (ZERO, ZERO), (ZERO, ZERO)))
    triangle_A = block_vector(
        0, ((ONE, ONE), (W, W), (W2, W2), (ZERO, ZERO))
    )
    triangle_B = block_vector(
        0, ((W, W), (ONE, ONE), (W2, W2), (ZERO, ZERO))
    )
    return {"secant": (secant_A, secant_B),
            "triangle": (triangle_A, triangle_B)}


def embed_selected_slices(A, B, C=None, D=None):
    """Embed q4 slices at ordered ports (0,1) into a six-port tensor."""
    slices = {(0, 0): A, (1, 1): B}
    if C is not None:
        slices[(0, 1)] = C
    if D is not None:
        slices[(1, 0)] = D
    f = [ZERO] * 64
    for endpoint, vector in slices.items():
        for y, value in zip(BITS4, vector):
            x = endpoint + y
            f[INDEX6[x]] = value
    return tuple(f)


def card(f, pair, sector, phase):
    """q_s=f(00)-s f(11) or f(01)-s f(10), in residual port order."""
    i, j = pair
    endpoints = ((0, 0), (1, 1)) if sector == 0 else ((0, 1), (1, 0))
    answer = []
    for y in BITS4:
        values = []
        for endpoint in endpoints:
            x = []
            y_it = iter(y)
            for k in range(6):
                if k == i:
                    x.append(endpoint[0])
                elif k == j:
                    x.append(endpoint[1])
                else:
                    x.append(next(y_it))
            values.append(f[INDEX6[tuple(x)]])
        answer.append(values[0] - phase * values[1])
    return tuple(answer)


def in_line(v, line_vector):
    if first_nonzero(v) is None:
        return True
    p = first_nonzero(line_vector)
    if p is None:
        raise AssertionError("zero target line")
    scalar = v[p] / line_vector[p]
    return v == vscale(scalar, line_vector)


ALL_SCALARS = object()


def allowed_scalars(q0, q1, atlas_vectors):
    """lambda with q0+lambda*q1 zero or on an atlas projective line."""
    answers = set()
    for v in atlas_vectors:
        q1_on_line = in_line(q1, v)
        if q1_on_line:
            if in_line(q0, v):
                return ALL_SCALARS
            continue
        p = first_nonzero(v)
        candidate = None
        for i in range(len(v)):
            denominator = q1[i] * v[p] - q1[p] * v[i]
            if denominator:
                numerator = q0[i] * v[p] - q0[p] * v[i]
                candidate = -numerator / denominator
                break
        assert candidate is not None
        if in_line(vadd(q0, vscale(candidate, q1)), v):
            answers.add(candidate)
    return answers


INCONSISTENT = object()


def canonical_system(rows):
    """RREF for affine equations a*x+b*y=c over Q(w)."""
    matrix = [list(row) for row in rows if row[0] or row[1] or row[2]]
    pivot_row = 0
    for column in (0, 1):
        pivot = next(
            (r for r in range(pivot_row, len(matrix)) if matrix[r][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        z = ONE / matrix[pivot_row][column]
        matrix[pivot_row] = [z * value for value in matrix[pivot_row]]
        for r in range(len(matrix)):
            if r == pivot_row:
                continue
            z = matrix[r][column]
            if z:
                matrix[r] = [
                    x - z * y for x, y in zip(matrix[r], matrix[pivot_row])
                ]
        pivot_row += 1
    for row in matrix:
        if not row[0] and not row[1] and row[2]:
            return INCONSISTENT
    nonzero = [tuple(row) for row in matrix if row[0] or row[1]]
    nonzero.sort(key=lambda row: 0 if row[0] else 1)
    return tuple(nonzero)


def intersect_systems(left, right):
    return canonical_system(left + right)


def prune_union(systems):
    unique = []
    for system in systems:
        if system is INCONSISTENT or system in unique:
            continue
        unique.append(system)
    keep = []
    for i, system in enumerate(unique):
        contained = False
        for j, other in enumerate(unique):
            if i == j:
                continue
            # system subset other iff adding other's equations changes nothing.
            if intersect_systems(system, other) == system:
                contained = True
                break
        if not contained:
            keep.append(system)
    return tuple(keep)


def affine_preimage_system(q0, q1, q2, target):
    """Equations for q0+c*q1+d*q2 to lie on span(target)."""
    p = first_nonzero(target)
    rows = []
    for i in range(len(target)):
        if i == p:
            continue
        a = q1[i] * target[p] - q1[p] * target[i]
        b = q2[i] * target[p] - q2[p] * target[i]
        constant = -(q0[i] * target[p] - q0[p] * target[i])
        rows.append((a, b, constant))
    return canonical_system(rows)


def preimage_union(q0, q1, q2, atlas_vectors):
    return prune_union(
        [affine_preimage_system(q0, q1, q2, target)
         for target in atlas_vectors]
    )


def system_rank(system):
    return len(system)


def point_of_system(system):
    """One deterministic point; for a line choose its free coordinate 0."""
    if len(system) == 0:
        return ZERO, ZERO
    if len(system) == 2:
        values = [ZERO, ZERO]
        for row in system:
            pivot = 0 if row[0] else 1
            values[pivot] = row[2]
        return tuple(values)
    row = system[0]
    if row[0]:
        return row[2], ZERO
    return ZERO, row[2]


def constant_switching_completion_audit(line_data):
    """Exhaust constant-projective opposite pencils, including both scalars."""
    odd_vectors = tuple(line_data[1][0])
    even_vectors = tuple(line_data[0][0])
    order = card_order()
    report = {}
    for name, (A, B) in switching_representatives().items():
        f0 = embed_selected_slices(A, B)
        final_components = []
        first_card_histogram = Counter()
        for direction_index, direction in enumerate(odd_vectors):
            f_c = embed_selected_slices(
                (ZERO,) * 16, (ZERO,) * 16, direction, (ZERO,) * 16
            )
            f_d = embed_selected_slices(
                (ZERO,) * 16, (ZERO,) * 16, (ZERO,) * 16, direction
            )
            current = ((),)
            for step, constraint in enumerate(order):
                pair, sector, phase = constraint
                targets = odd_vectors if sector == 1 else even_vectors
                allowed = preimage_union(
                    card(f0, pair, sector, phase),
                    card(f_c, pair, sector, phase),
                    card(f_d, pair, sector, phase),
                    targets,
                )
                intersections = []
                for old in current:
                    for new in allowed:
                        intersection = intersect_systems(old, new)
                        if intersection is not INCONSISTENT:
                            intersections.append(intersection)
                current = prune_union(intersections)
                if not current:
                    first_card_histogram[constraint] += 1
                    break
            final_components.extend((direction_index, system) for system in current)

        # Any component not equal to the isolated origin would be a completion.
        nonzero_components = []
        for direction_index, system in final_components:
            if system_rank(system) < 2:
                nonzero_components.append((direction_index, system))
                continue
            if point_of_system(system) != (ZERO, ZERO):
                nonzero_components.append((direction_index, system))
        assert not nonzero_components
        report[name] = {
            "directions": len(odd_vectors),
            "final_components": len(final_components),
            "first_failure_histogram": first_card_histogram,
        }

        # The identically zero opposite sector is c=d=0 and must itself fail.
        assert not all(
            atlas_contains(
                card(f0, pair, sector, phase),
                set(q4_atlas(sector)),
            )
            for pair in itertools.combinations(range(6), 2)
            for sector in (0, 1)
            for phase in MU3
        )
    return report


def card_order():
    """Deterministic exhaustive order, with the paper's early cards first."""
    preferred = [
        ((0, 2), 0, W), ((0, 3), 0, ONE), ((0, 3), 0, W),
        ((0, 2), 1, ONE), ((0, 3), 1, ONE),
        ((0, 3), 1, W), ((0, 2), 1, W),
    ]
    all_cards = [
        (pair, sector, phase)
        for pair in itertools.combinations(range(6), 2)
        if pair != (0, 1)
        for sector in (0, 1)
        for phase in MU3
    ]
    return tuple(preferred + [x for x in all_cards if x not in preferred])


def validate_switching_representatives(line_data):
    even_atlas, even_rich, _ = line_data[0]
    even_set = set(even_atlas)
    reps = switching_representatives()
    for name, (A, B) in reps.items():
        assert all(atlas_contains(vsub(A, vscale(s, B)), even_set) for s in MU3)
        line = line_key(vsub(A, B), vsub(A, vscale(W, B)))
        ids = even_rich[line]
        types = tuple(sorted(even_atlas[tuple(even_atlas)[i]][0] for i in ids))
        if name == "secant":
            assert types == (1, 1, 2)
        else:
            assert types == (2, 2, 2)


def nonconstant_switching_completion_audit(line_data):
    """Exhaust all nonconstant opposite pencils and their relative scale."""
    odd_atlas, odd_rich, _ = line_data[1]
    odd_vectors = tuple(odd_atlas)
    opposite = labeled_pencils(odd_atlas, odd_rich)
    assert len(opposite) == 202 * 6 == 1212
    order = card_order()
    report = {}
    for name, (A, B) in switching_representatives().items():
        f0 = embed_selected_slices(A, B)
        failure_histogram = Counter()
        total_scalar_candidates = 0
        completions = []
        for C, D in opposite:
            f1 = embed_selected_slices(
                (ZERO,) * 16, (ZERO,) * 16, C, D
            )
            possible = ALL_SCALARS
            failed_at = None
            for constraint in order:
                pair, sector, phase = constraint
                q0 = card(f0, pair, sector, phase)
                q1 = card(f1, pair, sector, phase)
                allowed = allowed_scalars(q0, q1, odd_vectors if (
                    # Total parity is even: residual parity equals sector.
                    sector == 1
                ) else tuple(q4_atlas(0)))
                if allowed is ALL_SCALARS:
                    continue
                if possible is ALL_SCALARS:
                    possible = set(allowed)
                else:
                    possible &= allowed
                if not possible:
                    failed_at = constraint
                    break
            if possible is ALL_SCALARS:
                raise AssertionError("unconstrained continuum completion")
            possible.discard(ZERO)  # lambda=0 is the separately checked zero sector.
            total_scalar_candidates += len(possible)
            if possible:
                for scalar in possible:
                    f = vadd(f0, vscale(scalar, f1))
                    assert all(
                        atlas_contains(
                            card(f, pair, sector, phase),
                            set(q4_atlas(sector)),
                        )
                        for pair in itertools.combinations(range(6), 2)
                        for sector in (0, 1)
                        for phase in MU3
                    )
                    completions.append((C, D, scalar))
            else:
                failure_histogram[failed_at] += 1
        assert not completions
        report[name] = {
            "opposite_pencils": len(opposite),
            "surviving_scalar_candidates": total_scalar_candidates,
            "first_failure_histogram": failure_histogram,
        }
    return report


def residual_index(full_bits, deleted_pair):
    return tuple(full_bits[k] for k in range(6) if k not in deleted_pair)


def sector_union_mask(support_mask, pair, sector):
    """Union of the two endpoint-slice supports, as a 16-bit mask."""
    i, j = pair
    answer = 0
    for word_index, x in enumerate(BITS6):
        if not (support_mask >> word_index) & 1:
            continue
        if (x[i] ^ x[j]) != sector:
            continue
        answer |= 1 << INDEX4[residual_index(x, pair)]
    return answer


def allowed_sector_unions():
    answer = set()
    for parity in (0, 1):
        blocks = BLOCKS4[parity]
        for count in (1, 2):
            for selected in itertools.combinations(blocks, count):
                mask = 0
                for block in selected:
                    for x in block:
                        mask |= 1 << INDEX4[x]
                answer.add(mask)
    assert len(answer) == 8 + 12
    return frozenset(answer)


ALLOWED_UNIONS = allowed_sector_unions()


def local_sector_choices(pair, sector):
    """All global 64-bit supports restricted to one pair-sector.

    For a chosen union U of one/two complement blocks, every word of U is
    placed in slice 0 only, slice 1 only, or both.  This is exactly the
    1+8*3^2+12*3^4 count, including the empty sector.
    """
    i, j = pair
    endpoints = ((0, 0), (1, 1)) if sector == 0 else ((0, 1), (1, 0))
    result = [0]
    for union in sorted(ALLOWED_UNIONS):
        residual_words = [BITS4[k] for k in range(16) if (union >> k) & 1]
        for placements in itertools.product((1, 2, 3), repeat=len(residual_words)):
            mask = 0
            for y, placement in zip(residual_words, placements):
                for endpoint_number, endpoint in enumerate(endpoints):
                    if not placement & (1 << endpoint_number):
                        continue
                    x = []
                    y_it = iter(y)
                    for k in range(6):
                        if k == i:
                            x.append(endpoint[0])
                        elif k == j:
                            x.append(endpoint[1])
                        else:
                            x.append(next(y_it))
                    mask |= 1 << INDEX6[tuple(x)]
            result.append(mask)
    assert len(result) == 1 + 8 * 3**2 + 12 * 3**4 == 1045
    return result


def even_partitions(items):
    """Set partitions into even blocks, each of size at least two."""
    items = tuple(items)
    if not items:
        yield ()
        return
    first = items[0]
    rest = items[1:]
    for size in range(2, len(items) + 1, 2):
        for companions in itertools.combinations(rest, size - 1):
            block = (first,) + companions
            remaining = tuple(x for x in rest if x not in companions)
            for tail in even_partitions(remaining):
                yield (block,) + tail


def affine_partition_supports():
    classes = Counter()
    supports = set()
    for partition in even_partitions(range(6)):
        direction = []
        for choice in itertools.product((0, 1), repeat=len(partition)):
            x = [0] * 6
            for use, block in zip(choice, partition):
                if use:
                    for i in block:
                        x[i] = 1
            direction.append(tuple(x))
        unseen = set(BITS6)
        while unseen:
            a = min(unseen)
            coset = {tuple(u ^ v for u, v in zip(a, d)) for d in direction}
            unseen -= coset
            mask = sum(1 << INDEX6[x] for x in coset)
            supports.add(mask)
            classes[tuple(sorted(map(len, partition), reverse=True))] += 1
    assert classes == Counter({(4, 2): 240, (2, 2, 2): 120, (6,): 32})
    assert len(supports) == 392
    return supports, classes


def permute_word(x, permutation):
    return tuple(x[permutation[i]] for i in range(6))


def nonproduct_orbit():
    decimal_words = (0, 3, 12, 51, 60, 63)
    # Decimal notation in the paper is standard six-bit big endian.
    base_words = {
        tuple((value >> (5 - i)) & 1 for i in range(6))
        for value in decimal_words
    }
    orbit = set()
    for permutation in itertools.permutations(range(6)):
        permuted = [permute_word(x, permutation) for x in base_words]
        for flip in BITS6:
            words = {
                tuple(a ^ b for a, b in zip(x, flip)) for x in permuted
            }
            orbit.add(sum(1 << INDEX6[x] for x in words))
    assert len(orbit) == 480
    return orbit


def support_audit():
    pair0 = (0, 1)
    equal = local_sector_choices(pair0, 0)
    mixed = local_sector_choices(pair0, 1)
    survivors = []
    pairs = tuple(itertools.combinations(range(6), 2))
    for a in equal:
        for b in mixed:
            support = a | b
            if all(
                sector_union_mask(support, pair, sector) in ALLOWED_UNIONS
                or sector_union_mask(support, pair, sector) == 0
                for pair in pairs[1:]
                for sector in (0, 1)
            ):
                survivors.append(support)
    assert len(survivors) == len(set(survivors)) == 873
    products, classes = affine_partition_supports()
    exceptional = set(survivors) - products - {0}
    assert len(set(survivors) & products) == 392
    assert len(exceptional) == 480
    orbit = nonproduct_orbit()
    assert exceptional == orbit
    return {
        "local_equal": len(equal),
        "local_mixed": len(mixed),
        "survivors": len(survivors),
        "product_classes": classes,
        "exceptional": len(exceptional),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-support", action="store_true",
        help=(
            "skip the 1045^2 all-pair support intersection; the final status "
            "is then PARTIAL PASS rather than PASS"
        ),
    )
    switching = parser.add_mutually_exclusive_group()
    switching.add_argument(
        "--switching",
        dest="switching",
        action="store_true",
        default=True,
        help=(
            "exhaust nonconstant and constant opposite pencils for the two "
            "switches (the default; retained for command compatibility)"
        ),
    )
    switching.add_argument(
        "--skip-switching",
        dest="switching",
        action="store_false",
        help=(
            "skip the switching-completion audit; the final status is then "
            "PARTIAL PASS rather than PASS"
        ),
    )
    args = parser.parse_args()
    started = time.monotonic()
    lines = q4_line_audit()
    print("q4 atlas: 66 points per parity, 132 total")
    print("rich lines per parity:", dict(lines[0][2]))
    print("labeled nonconstant pencils: 2 * 202 * 6 = 2424")
    validate_switching_representatives(lines)
    if not args.skip_support:
        result = support_audit()
        print("one-sector choices:", result["local_equal"])
        print("all-pair support survivors:", result["survivors"])
        print("product supports:", dict(result["product_classes"]))
        print("exceptional orbit:", result["exceptional"])
    if args.switching:
        report = nonconstant_switching_completion_audit(lines)
        for name, row in report.items():
            print(name, "nonconstant opposite pencils:", row["opposite_pencils"])
            print(name, "surviving relative scales:",
                  row["surviving_scalar_candidates"])
            print(name, "first-failure histogram:",
                  dict(row["first_failure_histogram"]))
        constant_report = constant_switching_completion_audit(lines)
        for name, row in constant_report.items():
            print(name, "constant directions:", row["directions"])
            print(name, "constant final affine components:",
                  row["final_components"])
    skipped = []
    if args.skip_support:
        skipped.append("all-pair support intersection")
    if not args.switching:
        skipped.append("switching-completion audit")
    if skipped:
        print("SKIPPED: " + "; ".join(skipped))
        print(f"PARTIAL PASS ({time.monotonic() - started:.3f}s)")
    else:
        print(f"PASS ({time.monotonic() - started:.3f}s)")


if __name__ == "__main__":
    main()
