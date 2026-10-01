"""Verify every genus-two extreme ray inside GL_2(F_3).

Run from any directory: python verify_bolza_witnesses.py
The identification of this group with the automorphism group of the Bolza
curve is the cited geometric input. All subgroup and index checks here
use exact integer arithmetic.
"""

import json
from pathlib import Path

import numpy as np

from genus2_compute import groups, enumerate_subgroups, profile, canonical


def main():
    data = Path(__file__).resolve().parent / "genus2_data"
    payload = json.loads((data / "bolza_ray_witnesses.json").read_text())
    assert payload["carrier"] == "GL2F3"
    group = next(group for group in groups() if group["name"] == "GL2F3")
    subgroups = set(enumerate_subgroups(group))
    assert len(group["mul"]) == 48
    assert len(subgroups) == 55
    expected = {
        tuple(map(int, vector))
        for vector in np.load(data / "extreme_rays_verified.npz")["indices"]
    }
    assert len(expected) == 127
    witnesses = payload["rays"]
    indices = [tuple(witness["indices"]) for witness in witnesses]
    assert len(indices) == len(set(indices)) == 127
    assert set(indices) == expected
    for witness in witnesses:
        ambient = witness["ambient_bitmask"]
        observation_subgroups = witness["subgroups_bitmasks"]
        assert ambient in subgroups
        assert len(observation_subgroups) == 4
        assert all(
            subgroup in subgroups and subgroup & ambient == subgroup
            for subgroup in observation_subgroups
        )
        assert profile(
            observation_subgroups, ambient.bit_count()
        ) == tuple(witness["indices"])
    assert len({canonical(vector) for vector in expected}) == 18
    print("PASS: GL2(F3) group axioms and all 55 subgroups regenerated.")
    print("PASS: all 127 extreme-ray indices have exact subgroup witnesses.")
    print("PASS: these witnesses realize all 18 variable-permutation orbits.")


if __name__ == "__main__":
    main()
