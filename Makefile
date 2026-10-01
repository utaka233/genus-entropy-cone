.PHONY: verify

verify:
	cd supplement && python3 verify_genus2.py
	cd supplement && python3 verify_bolza_witnesses.py
	cd supplement && python3 genus2_other_face_ranks.py
	cd supplement && python3 uniform_outer_verify.py
	cd supplement && python3 uniform_outer_dfz_verify.py
	cd supplement && python3 uniform_modulus_verify.py
	cd supplement && python3 group_structure_verify.py
	cd supplement && python3 verify_genus3_obstruction.py
	cd supplement && python3 group_ingleton_verify.py
	python3 revision/score_witness.py
