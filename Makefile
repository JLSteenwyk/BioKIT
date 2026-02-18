PYTHON := $(if $(wildcard venv/bin/python),venv/bin/python,python3)

## Run alignment-based functions
run.alignment_summary:
	python3 -m biokit-runner alignment_summary ./tests/sample_files/simple.fa 

## Install, develop, and testing make commands
dev.setup:
	python3 -m pip install --upgrade pip
	python3 -m pip install -r requirements-dev.txt
	python3 -m pip install -e .

install:
	# install so biokit command is available in terminal
	python3 setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python setup.py develop

test: test.unit test.integration

test.unit:
	python3 -m pytest -m "not integration"

test.integration:
	rm -rf output/
	mkdir output/
	python3 -m pytest --basetemp=output -m "integration"

test.fast:
	python -m pytest -m "not (integration or slow)"
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output -m "integration and not slow"

lint:
	python3 -m flake8 biokit tests

typecheck:
	python3 -m mypy biokit/helpers biokit/services/coding_sequences biokit/services/alignment biokit/services/text

security:
	python3 -m pip_audit

guard.private-api:
	bash scripts/check_private_api_usage.sh

quality: lint typecheck guard.private-api test.fast

bench.profile:
	$(PYTHON) scripts/benchmark_hotspots.py --python $(PYTHON) --iterations $${N:-5} --warmups $${WARMUPS:-1}

bench.profile.save:
	mkdir -p output
	$(PYTHON) scripts/benchmark_hotspots.py --python $(PYTHON) --iterations $${N:-5} --warmups $${WARMUPS:-1} --output output/bench_profile.tsv

precommit.install:
	python3 -m pre_commit install

precommit.run:
	python3 -m pre_commit run --all-files

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
