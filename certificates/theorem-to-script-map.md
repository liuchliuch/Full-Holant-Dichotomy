# Theorem-to-script map

The labels below are the stable LaTeX labels in the manuscript.
A script checks the specified finite claim, not every argument in the theorem
that uses it. Requirements and commands are in the [README](README.md).

| Finite claim | Program | Use in the paper |
| --- | --- | --- |
| Proposition A.1 (`prop:cert-odd-dih-q6-order-three`): order-three projective quaternary set, six-port exclusions, and all-pair support classification. | [`verify_odd_dihedral.py`](verify_odd_dihedral.py) | Theorem A.3 (`thm:cert-odd-dih-q6-interface`), then Theorem 5.14 (`thm:odd-dih-q6-base`). |
| Proposition A.4 (`prop:cert-high-charge-q6-joint-deck`): the implication for the displayed six-port loop and two-copy arrays. | [`verify_high_charge_q6.py`](verify_high_charge_q6.py) | Theorem A.6 (`thm:cert-high-charge-q6-interface`) and Lemma `lem:hc-q6-joint-rigidity`; used in Theorem 7.5 (`thm:hc-q6-complete`) and the height-three argument. |
| Proposition A.8 (`prop:cert-q8-collision-free-atlas`): nine matching families, seven induced by global perfect matchings and two of translation form. | [`verify_high_charge_q8_atlas.py`](verify_high_charge_q8_atlas.py) | Theorem A.9 (`thm:cert-q8-collision-free-interface`), then Lemma 8.9 (`lem:q8-collision-free-atlas`) in the height-four argument. |
| Proposition A.10 (`prop:cert-q8-canonical-support`): the 14-plane/56-nonplane split, nonplane graph connectivity and odd cycle, and plane--port-pair incidence connectivity. | [`verify_high_charge_q8_atlas.py`](verify_high_charge_q8_atlas.py) | Supplementary checks of the printed proof; Theorem A.11 (`thm:cert-q8-canonical-support-interface`) and Lemma 8.10 (`lem:q8-canonical-support`) use its conclusion. |
| Nine 560-by-98 coefficient/scaling systems of rank 97 and the 672-by-14 multiplicative relation lattice with Smith data `1^9, 2^4, 0`. | [`verify_high_charge_q8_normalized.py`](verify_high_charge_q8_normalized.py) | Supplementary calculation; no separate proof dependency. |

In Proposition A.1(ii), nonconstancy is projective. The contraction hypothesis
ranges over every required port pair, xor sector, and phase. Support
classification alone does not establish coefficient factorization; that
argument remains in the text.

The six-port wrapper invokes two internal stages:
`verify_high_charge_q6_orbits.cpp` enumerates 14,348,907 signed branches,
1,868 linear spaces, and 24 orbits. `verify_high_charge_q6_algebra.py`
checks agreement of the orbit sets and makes 61 exact real-polynomial queries
(51 sparse charts and 10 dense implications) for the 23 nonzero families.
Its polynomial conditions are necessary for the quaternary matching-product
class in the paper. The UNSAT conclusions rely on Z3 4.14.0.

The normalized eight-port calculation assumes (HN.5) at eight ports: every
pair contraction is a nonzero scalar multiple of an X-matching product.
Its lattice classification concerns nonzero solutions of the generated
multiplicative relations. It does not establish the preceding normalization
or the subsequent reduction adjoining equality.
