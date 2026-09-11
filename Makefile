.PHONY: install verify verify-all experiments docs clean

install:
	python -m pip install -r requirements.txt

verify:
	python code/run_verifications.py --rational-only

verify-all:
	python code/run_verifications.py

experiments:
	python code/run_verifications.py --experiments

docs:
	python build_supplement.py

clean:
	rm -f supplementary_material_gamma.aux supplementary_material_gamma.log supplementary_material_gamma.out supplementary_material_gamma.toc
