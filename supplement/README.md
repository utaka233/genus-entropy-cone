# Exact verification suite

Run **make verify** from the repository root. The ten entries are listed
in the root Makefile. NumPy, SciPy and a C++17 compiler named c++ are required.

Coordinates are bitmasks 1,...,15:
1,2,12,3,13,23,123,4,14,24,124,34,134,234,1234.
An integer index vector d represents the entropy vector log d.

Input certificates are profiles_canonical.pkl, extreme_rays_verified.npz,
extreme_rays_certificates.json, lp_basis_certificates.pkl,
ingleton_genus_comparison.json, bolza_ray_witnesses.json,
uniform_outer_certificates.json and uniform_outer_dfz_certificate.json.
The genus-two data files are under genus2_data/. Keep their directory layout.

The numerical LP/QR calculations select candidate witnesses; membership,
exposure, rank and slice claims are checked with integer and rational
arithmetic. Other symmetric-face formal rays are derived from the already
verified integer ray indices. polynomial_tools.py contains only the exact
arithmetic helpers used by this rank check.

The scripts regenerate result JSON files. No result log is used as a
substitute for running a check. The genus-three and larger-group searches
construct their groups internally. revision/score_witness.py contains its
rational probabilities directly.

The published curve-action classification, the order-less-than-120 exclusion
and the best published score bound remain external inputs. Some older repair
identities verify valid weaker bounds; they are retained for the identities,
not presented as the manuscript's strongest current error bounds.
