.PHONY: install verify verify-all test experiments data-check docs integrity clean

install:
	python -m pip install -r requirements.txt

verify:
	python code/run_verifications.py --rational-only

verify-all:
	python code/run_verifications.py

test:
	python -m unittest discover -s tests -v

experiments:
	python code/run_verifications.py --experiments

data-check:
	python code/verify_experiment_data.py

docs:
	python build_supplement.py

integrity:
	python code/check_integrity.py

clean:
	python -c "import pathlib,shutil; [shutil.rmtree(p) for p in pathlib.Path('.').rglob('__pycache__')]"
