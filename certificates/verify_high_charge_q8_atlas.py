#!/usr/bin/env python3
"""Exact replay of the arity-eight high-charge matching/support atlas.

The calculation is purely finite and uses only Python's standard library.
It reconstructs every compatible derivative-matching deck after fixing one
matching by symmetry, separates the overlapping-factor decks from the two
Reed--Muller translation decks, and verifies the connected odd-cycle
elimination of all 56 nonplane four-subsets together with the complementary
plane--crossing-pair incidence connectivity.
"""

from __future__ import annotations

from collections import Counter, deque
from itertools import combinations, permutations

if not __debug__:
    raise SystemExit("exact verifier requires assertions; do not use python -O")


PORTS = tuple(range(8))
PAIRS = tuple(combinations(PORTS, 2))


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


OPTIONS = {
    deleted: tuple(
        perfect_matchings(tuple(v for v in PORTS if v not in deleted))
    )
    for deleted in PAIRS
}
ALL_MATCHINGS = tuple(perfect_matchings(PORTS))


def induced_matching(
    matching: tuple[tuple[int, int], ...],
    second_deleted: tuple[int, int],
) -> tuple[tuple[int, int], ...]:
    """Matching left after one further derivative contraction."""

    if second_deleted in matching:
        return tuple(edge for edge in matching if edge != second_deleted)
    a, b = second_deleted
    edge_a = next(edge for edge in matching if a in edge)
    edge_b = next(edge for edge in matching if b in edge)
    mate_a = edge_a[0] if edge_a[1] == a else edge_a[1]
    mate_b = edge_b[0] if edge_b[1] == b else edge_b[1]
    untouched = next(edge for edge in matching if edge not in (edge_a, edge_b))
    return tuple(sorted((tuple(sorted((mate_a, mate_b))), untouched)))


def compatible(
    deleted: tuple[int, int],
    matching: tuple[tuple[int, int], ...],
    chosen: dict[tuple[int, int], tuple[tuple[int, int], ...]],
) -> bool:
    for previous_deleted, previous_matching in chosen.items():
        if not disjoint(deleted, previous_deleted):
            continue
        if induced_matching(matching, previous_deleted) != induced_matching(
            previous_matching, deleted
        ):
            return False
    return True


def enumerate_decks():
    chosen = {(0, 1): ((2, 3), (4, 5), (6, 7))}
    solutions = []

    def search() -> None:
        if len(chosen) == len(PAIRS):
            solutions.append(dict(chosen))
            return
        best_pair = None
        best_options = None
        for deleted in PAIRS:
            if deleted in chosen:
                continue
            options = [
                matching
                for matching in OPTIONS[deleted]
                if compatible(deleted, matching, chosen)
            ]
            if not options:
                return
            if best_options is None or len(options) < len(best_options):
                best_pair, best_options = deleted, options
        assert best_pair is not None and best_options is not None
        for matching in best_options:
            chosen[best_pair] = matching
            search()
            del chosen[best_pair]

    search()
    return solutions


def overlap_count(deck) -> int:
    return sum(
        1
        for left, right in combinations(PAIRS, 2)
        if not disjoint(left, right) and set(deck[left]) & set(deck[right])
    )


def global_matching_card(matching, deleted):
    """Three-edge card induced by one global four-edge matching."""

    if deleted in matching:
        return tuple(edge for edge in matching if edge != deleted)
    a, b = deleted
    edge_a = next(edge for edge in matching if a in edge)
    edge_b = next(edge for edge in matching if b in edge)
    mate_a = edge_a[0] if edge_a[1] == a else edge_a[1]
    mate_b = edge_b[0] if edge_b[1] == b else edge_b[1]
    untouched = tuple(edge for edge in matching if edge not in (edge_a, edge_b))
    return tuple(sorted((tuple(sorted((mate_a, mate_b))), *untouched)))


def global_matching_ancestors(deck):
    return [
        matching
        for matching in ALL_MATCHINGS
        if all(global_matching_card(matching, edge) == deck[edge] for edge in PAIRS)
    ]


def translation_matching(deleted: tuple[int, int]):
    step = deleted[0] ^ deleted[1]
    seen = set(deleted)
    matching = []
    for vertex in PORTS:
        mate = vertex ^ step
        if vertex in seen or mate in seen:
            continue
        seen.update((vertex, mate))
        matching.append(tuple(sorted((vertex, mate))))
    return tuple(sorted(matching))


def is_translation_atlas_up_to_relabelling(deck) -> bool:
    for permutation in permutations(PORTS):
        valid = True
        for deleted, matching in deck.items():
            image_deleted = tuple(sorted(permutation[v] for v in deleted))
            image_matching = tuple(
                sorted(tuple(sorted(permutation[v] for v in edge)) for edge in matching)
            )
            if image_matching != translation_matching(image_deleted):
                valid = False
                break
        if valid:
            return True
    return False


def verify_nonplane_elimination() -> None:
    four_sets = tuple(frozenset(block) for block in combinations(PORTS, 4))
    planes = {block for block in four_sets if _xor(block) == 0}
    nonplanes = set(four_sets) - planes
    assert len(planes) == 14 and len(nonplanes) == 56

    graph = {block: set() for block in nonplanes}
    for block in nonplanes:
        sigma = _xor(block)
        for p in block:
            for q in set(PORTS) - set(block):
                step = p ^ q
                if sigma in (0, step):
                    continue
                swapped = frozenset((set(block) - {p}) | {q})
                if swapped in nonplanes:
                    graph[block].add(swapped)

    start = next(iter(nonplanes))
    reached = {start}
    queue = deque([start])
    while queue:
        block = queue.popleft()
        for neighbor in graph[block]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    assert reached == nonplanes

    triangle = tuple(
        frozenset(block)
        for block in ((0, 1, 2, 4), (0, 1, 2, 5), (0, 2, 4, 5))
    )
    assert all(triangle[(i + 1) % 3] in graph[triangle[i]] for i in range(3))

    plane_list = tuple(sorted(planes, key=lambda block: tuple(sorted(block))))
    nodes = [*(('plane', i) for i in range(len(plane_list))),
             *(('pair', i) for i in range(len(PAIRS)))]
    incidence = {node: set() for node in nodes}
    for plane_index, plane in enumerate(plane_list):
        for pair_index, pair in enumerate(PAIRS):
            if sum(vertex in plane for vertex in pair) == 1:
                plane_node = ('plane', plane_index)
                pair_node = ('pair', pair_index)
                incidence[plane_node].add(pair_node)
                incidence[pair_node].add(plane_node)
    reached = {nodes[0]}
    queue = deque([nodes[0]])
    while queue:
        node = queue.popleft()
        for neighbor in incidence[node]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    assert reached == set(nodes)


def _xor(values) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def main() -> None:
    decks = enumerate_decks()
    histogram = Counter(overlap_count(deck) for deck in decks)
    assert len(decks) == 9
    assert histogram == Counter({168: 7, 0: 2})

    overlapping = [deck for deck in decks if overlap_count(deck) != 0]
    ancestors = [global_matching_ancestors(deck) for deck in overlapping]
    assert all(len(matches) == 1 for matches in ancestors)

    collision_free = [deck for deck in decks if overlap_count(deck) == 0]
    assert len(collision_free) == 2
    assert all(is_translation_atlas_up_to_relabelling(deck) for deck in collision_free)

    verify_nonplane_elimination()
    print(
        "compatible derivative decks: 9 = 7 unique global-matching ancestors "
        "+ 2 translation atlases"
    )
    print(
        "canonical support atlas: 14 affine planes, 56 nonplanes eliminated; "
        "plane-pair incidence connected"
    )
    print("PASS")


if __name__ == "__main__":
    main()
