# Verification record

Prepared on 2026-10-01. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8; the exhaustive group search also uses a C++17 compiler.

## Explanatory notebook

All ten Python code cells of `demo.ipynb` were executed sequentially in a
fresh Python namespace. All assertions passed, and the three saved figures
were inspected. Outputs in the notebook come from this execution.

The preparation environment did not have Jupyter/IPython installed, so this
was a direct execution of the cell sources with the stated Python libraries,
not an execution through a Jupyter kernel. Matplotlib's Agg backend was used
to capture the figures as standard notebook image outputs. Starting
JupyterLab and invoking its Run All command have not been tested here.

The notebook uses only small examples and embedded teaching data. Its
floating-point calculations explain the results; they do not replace the
exact certificate checks below.

## Exact computational supplement

`make verify` completed successfully with the reduced package, including all
ten entry points:

1. `verify_genus2.py`: exact membership and extremality certificates for
   the 127 genus-two rays.
2. `verify_bolza_witnesses.py`: the recorded subgroup witnesses.
3. `genus2_other_face_ranks.py`: the stated invariant face ranks.
4. `uniform_outer_verify.py`: rational repair and Zhang--Yeung identities.
5. `uniform_outer_dfz_verify.py`: the additional rational DFZ identity.
6. `uniform_modulus_verify.py`: the four-atom bounds and the stated
   Riemann--Hurwitz arithmetic exclusions.
7. `group_structure_verify.py`: the order-24 carrier presentation and
   its semidirect-product structure.
8. `verify_genus3_obstruction.py`: the exact subgroup obstruction.
9. `group_ingleton_verify.py`: 259,532,100 PSL2(F7) and 149,964,516 S5
   subgroup-multiset cases. The former has no offenders; the latter has
   60, all with ratio 15/16 and trivial fourfold intersection.
10. `revision/score_witness.py`: forty rational atoms and rational interval
    certificates for the improved score and error lower bounds.

Some repair scripts report the weaker, valid bounds 1/60 or 7/440 and the
earlier four-atom lower bound. The notebook instead displays the stronger
current interval, approximately [0.007895330349, 0.015789473684]. Its lower
endpoint is certified by `score_witness.py`; its upper endpoint 3/190 uses
the published score theorem cited in the manuscript.

## External inputs and scope

The curve-action classification, exclusion of violating groups of order
less than 120, and the published global score bound are mathematical
inputs from the cited literature. The package verifies the finite
computations used with these inputs; it does not reproduce their proofs or
the complete SmallGroups census. It does not compute the full entropy cone
or a general convergence modulus. No TeX compilation is needed to run it.
