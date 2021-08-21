## Run alignment-based functions
run.alignment_summary:
	python3 -m biokit-runner alignment_summary ./tests/sample_files/simple.fa 

## Install, develop, and testing make commands
install:
	# install so biokit command is available in terminal
	python setup.py install

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
	rm test.fa test.occupancy test.partition

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
