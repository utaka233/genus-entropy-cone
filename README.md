# A Genus Filtration of the Entropy Cone

Takayuki Uchiba. Computational supplement, prepared version 1.0.0.

Start with **demo.ipynb** for a self-contained walkthrough of the paper:
coset-valued random variables, the genus-two rays and symmetric slice,
an explicit S5 Ingleton violation, and the normalized genus-error bounds.
Its saved outputs can be read without rerunning the notebook.

## Run the demo

The demo needs Python, NumPy, SciPy and Matplotlib. JupyterLab opens it
interactively. No C++ compiler or TeX installation is needed for the demo.

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-demo.txt
jupyter lab demo.ipynb
```

Use Run All to recompute the notebook. The examples are small; the exhaustive
subgroup and certificate checks are separate commands below.

## Verify the exact computations

```sh
python3 -m pip install -r requirements.txt
make verify
```

The full suite needs a C++17 compiler available as c++ and GNU Make, in
addition to Python/NumPy/SciPy. It builds its C++ search in a temporary
directory and writes JSON results under supplement/ and revision/.
The computational environment and completed checks are in VERIFICATION.md.

## Contents

- demo.ipynb: explanatory computations with saved outputs.
- supplement/: the verifier entries and their shared exact-arithmetic modules.
- supplement/genus2_data/: six input data/certificate files.
- supplement/uniform_outer_certificates.json and
  uniform_outer_dfz_certificate.json: rational repair certificates.
- CITATION.cff: author, software title and prepared version metadata.

Exploratory full-facet searches, candidate-generation scripts, redundant
intermediate datasets and historical logs are omitted. Verification logic
is retained. The small polynomial helper module was separated from an
unused exploratory driver. The orbit data needed for teaching are embedded
in the notebook.

## Mathematical scope

The genus-two cone has 127 extreme rays in 18 permutation orbits; the
complete V-description and the stated invariant/symmetric facets are
verified. A complete 15-dimensional H-description is not claimed.

The genus-three obstruction and PSL2(F7)/S5 subgroup checks use integer
arithmetic. The repair identities and score intervals use rational
arithmetic. The complete curve-action classification and exclusion of
violating groups of order less than 120 are published external inputs.
This package does not rerun the SmallGroups census.

The notebook distinguishes the raw forty-atom score from the improved
limiting score. The displayed interval for E4 is a proved bound, not its
exact value. The upper endpoint 3/190 uses a published score theorem; this
suite does not re-prove that external theorem. A general computable
convergence modulus is not supplied.

The computations and manuscript were developed with AI assistance. The
programs perform exact computational verification, not proof-assistant
formalization.

## Citation and license

CITATION.cff supplies the citation metadata. Cite the published release
URL or its version-specific archive DOI. Select a software license when
creating the GitHub repository; place the selected LICENSE at the root.
