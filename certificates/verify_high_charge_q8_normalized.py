#!/usr/bin/env python3
"""Standard-library exact certificate for the normalized high-charge q8 base.

This file is deliberately self-contained.  It generates all matchings,
derivative decks, ancestor equations, and weighted chart relations from their
definitions.  It uses no serialized survivor table, random choice, floating
point arithmetic, checksum, or external CAS.

Scope.  The 560-by-98 systems assume that an earlier *actual* port-extension
argument has already put every nonzero derivative in the form

    partial_e F = lambda_e X_{M_e},  lambda_e != 0.

The certificate proves the finite normalized conclusion: the seven overlap
decks reconstruct one global matching product and the two collision-free
decks reconstruct the punctured RM(1,3) tensor. The final lattice check
classifies the nonzero solutions of the generated multiplicative rank-one
relations as a common scalar times one of sixteen sign characters.
It does not prove the preceding physical
normalization, factor extraction, endpoint repair, or whole-language Bell
completion.
"""

from __future__ import annotations

from functools import reduce
from itertools import combinations, permutations, product
from operator import xor
from time import perf_counter

if not __debug__:
    raise SystemExit("exact verifier requires assertions; do not use python -O")


PORTS = tuple(range(8))
EDGES = tuple(combinations(PORTS, 2))
MODULUS = 1_000_000_007  # prime


def disjoint(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return not (set(left) & set(right))


def perfect_matchings(vertices: tuple[int, ...]):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for index in range(1, len(vertices)):
        mate = vertices[index]
        rest = vertices[1:index] + vertices[index + 1 :]
        for tail in perfect_matchings(rest):
            yield tuple(sorted(((first, mate),) + tail))


MATCHING_OPTIONS = {
    deleted: tuple(
        perfect_matchings(tuple(v for v in PORTS if v not in deleted))
    )
    for deleted in EDGES
}
ALL_MATCHINGS = tuple(perfect_matchings(PORTS))


def induced_matching(
    matching: tuple[tuple[int, int], ...],
    contracted: tuple[int, int],
) -> tuple[tuple[int, int], ...]:
    """Matching left after one further X-contraction."""

    if contracted in matching:
        return tuple(edge for edge in matching if edge != contracted)
    a, b = contracted
    edge_a = next(edge for edge in matching if a in edge)
    edge_b = next(edge for edge in matching if b in edge)
    mate_a = edge_a[0] if edge_a[1] == a else edge_a[1]
    mate_b = edge_b[0] if edge_b[1] == b else edge_b[1]
    untouched = next(edge for edge in matching if edge not in (edge_a, edge_b))
    return tuple(sorted((tuple(sorted((mate_a, mate_b))), untouched)))


def enumerate_decks():
    chosen = {(0, 1): ((2, 3), (4, 5), (6, 7))}
    solutions = []
    nodes = 0

    def compatible(deleted, matching) -> bool:
        return all(
            not disjoint(deleted, previous_deleted)
            or induced_matching(matching, previous_deleted)
            == induced_matching(previous_matching, deleted)
            for previous_deleted, previous_matching in chosen.items()
        )

    def search() -> None:
        nonlocal nodes
        nodes += 1
        if len(chosen) == len(EDGES):
            solutions.append(dict(chosen))
            return
        best_edge = None
        best_choices = None
        for edge in EDGES:
            if edge in chosen:
                continue
            choices = [
                matching
                for matching in MATCHING_OPTIONS[edge]
                if compatible(edge, matching)
            ]
            if not choices:
                return
            if best_choices is None or len(choices) < len(best_choices):
                best_edge, best_choices = edge, choices
        assert best_edge is not None and best_choices is not None
        for matching in best_choices:
            chosen[best_edge] = matching
            search()
            del chosen[best_edge]

    search()
    return nodes, solutions


def has_overlap(deck) -> bool:
    return any(
        not disjoint(left, right) and set(deck[left]) & set(deck[right])
        for left, right in combinations(EDGES, 2)
    )


def translation_matching(edge: tuple[int, int]):
    step = edge[0] ^ edge[1]
    seen = set(edge)
    answer = []
    for vertex in PORTS:
        mate = vertex ^ step
        if vertex in seen or mate in seen:
            continue
        seen.update((vertex, mate))
        answer.append(tuple(sorted((vertex, mate))))
    return tuple(sorted(answer))


def translation_relabelling(deck):
    for permutation in permutations(PORTS):
        valid = True
        for deleted, matching in deck.items():
            image_deleted = tuple(sorted(permutation[v] for v in deleted))
            image_matching = tuple(
                sorted(
                    tuple(sorted(permutation[v] for v in edge))
                    for edge in matching
                )
            )
            if image_matching != translation_matching(image_deleted):
                valid = False
                break
        if valid:
            return permutation
    return None


WEIGHT_FOUR = tuple(
    sum(1 << vertex for vertex in subset)
    for subset in combinations(PORTS, 4)
)
WEIGHT_FOUR_INDEX = {mask: index for index, mask in enumerate(WEIGHT_FOUR)}


def ancestor_system(deck) -> list[list[int]]:
    """The 560 equations in 70 ancestor and 28 card-scalar variables."""

    rows = []
    for card_index, deleted in enumerate(EDGES):
        i, j = deleted
        residual = tuple(v for v in PORTS if v not in deleted)
        matching = deck[deleted]
        for triple in combinations(residual, 3):
            mask = sum(1 << v for v in triple)
            row = [0] * 98
            row[WEIGHT_FOUR_INDEX[mask | (1 << i)]] = 1
            row[WEIGHT_FOUR_INDEX[mask | (1 << j)]] = 1
            if all(sum(v in triple for v in edge) == 1 for edge in matching):
                row[70 + card_index] = -1
            rows.append(row)
    assert len(rows) == 560 and all(len(row) == 98 for row in rows)
    return rows


def sparse_rank_mod(rows: list[list[int]], modulus: int) -> int:
    """Exact rank over F_modulus by sparse Gaussian elimination."""

    basis: dict[int, dict[int, int]] = {}
    for source in rows:
        row = {
            column: value % modulus
            for column, value in enumerate(source)
            if value % modulus
        }
        while row:
            pivot = min(row)
            previous = basis.get(pivot)
            if previous is None:
                inverse = pow(row[pivot], -1, modulus)
                row = {
                    column: (value * inverse) % modulus
                    for column, value in row.items()
                    if (value * inverse) % modulus
                }
                basis[pivot] = row
                break
            factor = row[pivot]
            for column, value in previous.items():
                updated = (row.get(column, 0) - factor * value) % modulus
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
    return len(basis)


def transversal_support(matching) -> frozenset[int]:
    return frozenset(
        mask
        for mask in WEIGHT_FOUR
        if all(((mask >> a) & 1) ^ ((mask >> b) & 1) for a, b in matching)
    )


def candidate_vector(deck, support: frozenset[int]):
    """Return the exact 98-vector if support gives the prescribed X deck."""

    ancestor = [int(mask in support) for mask in WEIGHT_FOUR]
    scalars = []
    for deleted in EDGES:
        i, j = deleted
        residual = tuple(v for v in PORTS if v not in deleted)
        matching = deck[deleted]
        live_values = []
        for triple in combinations(residual, 3):
            mask = sum(1 << v for v in triple)
            value = (
                ancestor[WEIGHT_FOUR_INDEX[mask | (1 << i)]]
                + ancestor[WEIGHT_FOUR_INDEX[mask | (1 << j)]]
            )
            is_transversal = all(
                sum(v in triple for v in edge) == 1 for edge in matching
            )
            if is_transversal:
                live_values.append(value)
            elif value != 0:
                return None
        if not live_values or len(set(live_values)) != 1 or live_values[0] == 0:
            return None
        scalars.append(live_values[0])
    vector = ancestor + scalars
    assert len(vector) == 98
    return vector


def dot(left: list[int], right: list[int]) -> int:
    return sum(a * b for a, b in zip(left, right))


def xor_of_image(mask: int, permutation: tuple[int, ...]) -> int:
    return reduce(
        xor,
        (permutation[v] for v in PORTS if (mask >> v) & 1),
        0,
    )


def verify_ancestor_reconstruction(decks) -> tuple[int, int]:
    matching_decks = 0
    reed_muller_decks = 0
    for deck in decks:
        rows = ancestor_system(deck)
        # A nonzero exact kernel vector below proves rank_Q <= 97.  Rank 97
        # modulo one prime exhibits a 97-minor not divisible by that prime,
        # hence a nonzero integer 97-minor and rank_Q >= 97.  Therefore this
        # is an exact characteristic-zero rank proof, not a probabilistic one.
        assert sparse_rank_mod(rows, MODULUS) == 97

        if has_overlap(deck):
            witnesses = []
            for matching in ALL_MATCHINGS:
                support = transversal_support(matching)
                vector = candidate_vector(deck, support)
                if vector is not None:
                    witnesses.append((matching, support, vector))
            assert len(witnesses) == 1
            _, support, vector = witnesses[0]
            assert len(support) == 16
            matching_decks += 1
        else:
            permutation = translation_relabelling(deck)
            assert permutation is not None
            support = frozenset(
                mask
                for mask in WEIGHT_FOUR
                if xor_of_image(mask, permutation) == 0
            )
            assert len(support) == 14
            vector = candidate_vector(deck, support)
            assert vector is not None
            # No collision-free deck is simultaneously a matching cube.
            assert not any(
                candidate_vector(deck, transversal_support(matching)) is not None
                for matching in ALL_MATCHINGS
            )
            reed_muller_decks += 1

        assert all(dot(row, vector) == 0 for row in rows)
        assert all(value != 0 for value in vector[70:])

    assert (matching_decks, reed_muller_decks) == (7, 2)
    return matching_decks, reed_muller_decks


POINTS = tuple(product((0, 1), repeat=3))


def weighted_atlas_rows():
    """Generate all multiplicative rank-one relations on the 14 planes."""

    logical = []
    logical_index = {}
    for a in (0, 1):
        for u in POINTS:
            mask = sum(
                (
                    (a + sum(u[k] * point[k] for k in range(3))) & 1
                )
                << index
                for index, point in enumerate(POINTS)
            )
            if mask.bit_count() == 4:
                logical_index[(a, u)] = len(logical)
                logical.append((a, u))
    assert len(logical) == 14

    rows = []
    for p, q in EDGES:
        h = tuple(POINTS[p][k] ^ POINTS[q][k] for k in range(3))
        directions = [
            u
            for u in POINTS
            if sum(u[k] * h[k] for k in range(3)) & 1
        ]
        residual = [r for r in PORTS if r not in (p, q)]
        pairs = []
        seen = set()
        for r in residual:
            if r in seen:
                continue
            mate_point = tuple(POINTS[r][k] ^ h[k] for k in range(3))
            s = POINTS.index(mate_point)
            seen.update((r, s))
            pairs.append((r, s))
        assert len(pairs) == 3

        table = {}
        for a in (0, 1):
            for u in directions:
                y = tuple(
                    (
                        a
                        + sum(u[k] * POINTS[r][k] for k in range(3))
                    )
                    & 1
                    for r, _ in pairs
                )
                table[y] = logical_index[(a, u)]
        assert len(table) == 8

        for base in product((0, 1), repeat=3):
            for i, j in combinations(range(3), 2):
                corners = []
                for x, y in ((0, 0), (1, 1), (1, 0), (0, 1)):
                    point = list(base)
                    point[i], point[j] = x, y
                    corners.append(table[tuple(point)])
                row = [0] * 14
                row[corners[0]] += 1
                row[corners[1]] += 1
                row[corners[2]] -= 1
                row[corners[3]] -= 1
                rows.append(row)
    assert len(rows) == 672 and all(len(row) == 14 for row in rows)
    return logical, rows


def rank_mod_two(rows: list[list[int]], number_of_columns: int) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        mask = sum((value & 1) << index for index, value in enumerate(row))
        while mask:
            pivot = mask.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = mask
                break
            mask ^= basis[pivot]
    return len(basis)


def independent_rows_mod(
    rows: list[list[int]], columns: tuple[int, ...], modulus: int
) -> list[int]:
    basis: dict[int, dict[int, int]] = {}
    selected = []
    for row_index, source in enumerate(rows):
        row = {
            local: source[column] % modulus
            for local, column in enumerate(columns)
            if source[column] % modulus
        }
        while row:
            pivot = min(row)
            previous = basis.get(pivot)
            if previous is None:
                inverse = pow(row[pivot], -1, modulus)
                row = {
                    column: (value * inverse) % modulus
                    for column, value in row.items()
                    if (value * inverse) % modulus
                }
                basis[pivot] = row
                selected.append(row_index)
                break
            factor = row[pivot]
            for column, value in previous.items():
                updated = (row.get(column, 0) - factor * value) % modulus
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
        if len(selected) == len(columns):
            return selected
    return selected


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant with row pivoting."""

    size = len(matrix)
    assert size and all(len(row) == size for row in matrix)
    work = [row[:] for row in matrix]
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, size)
                if work[row][pivot_index] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row],
                work[pivot_index],
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % denominator == 0
                work[row][column] = numerator // denominator
            work[row][pivot_index] = 0
        denominator = pivot
    return sign * work[-1][-1]


def verify_weighted_character_group() -> None:
    logical, rows = weighted_atlas_rows()
    assert all(sum(row) == 0 for row in rows)

    # The row lattice L has rational rank 13: the all-one vector is in its
    # right kernel, while the deterministically selected 13-by-13 minor below
    # has determinant 16.  The torsion order of Z^14/L is the gcd of all
    # maximal minors, hence divides 16.  Modulo 2 the row rank is 9, so in
    # Smith language exactly four of the thirteen nonzero invariant factors
    # are even.  Their product is therefore divisible by 2^4.  The two facts
    # force the invariant factors to be 1^9, 2^4, followed by one free zero.
    columns = tuple(range(13))
    selected = independent_rows_mod(rows, columns, MODULUS)
    assert len(selected) == 13
    minor = [[rows[row][column] for column in columns] for row in selected]
    determinant = bareiss_determinant(minor)
    assert abs(determinant) == 16
    assert rank_mod_two(rows, 14) == 9

    constant = [1] * 14
    characters = [
        [a for a, _ in logical],
        *[[u[coordinate] for _, u in logical] for coordinate in range(3)],
    ]
    expected_kernel = [constant, *characters]
    assert rank_mod_two(expected_kernel, 14) == 5
    assert all(
        sum(coefficient * exponent for coefficient, exponent in zip(row, char))
        % 2
        == 0
        for row in rows
        for char in expected_kernel
    )

    # Since dim ker(A mod 2)=14-9=5, the five displayed vectors are the full
    # sign kernel.  The free factor is the common scalar (all relations have
    # degree zero); modulo that scalar the remaining sixteen solutions are
    # exactly (-1)^(alpha*a+beta.u).  Together with the invariant-factor
    # argument above this also excludes higher roots of unity and additional
    # continuous weighted deformations over C*.


def main() -> None:
    started = perf_counter()
    nodes, decks = enumerate_decks()
    assert (nodes, len(decks)) == (236, 9)

    overlap = [deck for deck in decks if has_overlap(deck)]
    collision_free = [deck for deck in decks if not has_overlap(deck)]
    assert (len(overlap), len(collision_free)) == (7, 2)
    translation_maps = [translation_relabelling(deck) for deck in collision_free]
    assert translation_maps == [
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 1, 2, 3, 4, 5, 7, 6),
    ]
    deck_time = perf_counter()

    matching_count, reed_muller_count = verify_ancestor_reconstruction(decks)
    ancestor_time = perf_counter()
    verify_weighted_character_group()
    finished = perf_counter()

    print("high-charge q8 normalized certificate: PASS")
    print(f"deck recursion: {nodes} nodes -> 9 decks = 7 overlap + 2 translation")
    print(
        "ancestor systems: nine 560x98 systems, exact rank 97; "
        f"{matching_count} unique 16-transversal + "
        f"{reed_muller_count} unique 14-plane"
    )
    print(
        "weighted atlas: 672x14 relations; equivalent Smith data "
        "1^9,2^4,0; characters C* x (mu_2)^4"
    )
    print(
        "timing seconds: "
        f"decks={deck_time - started:.3f}, "
        f"ancestors={ancestor_time - deck_time:.3f}, "
        f"weighted={finished - ancestor_time:.3f}, "
        f"total={finished - started:.3f}"
    )


if __name__ == "__main__":
    main()
