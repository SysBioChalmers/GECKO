.PHONY: clean clean-test clean-pyc clean-build docs help copy-models protein-properties
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts


clean-build: ## remove build artifacts
	rm -rf build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -type f -name '*.py[co]' -delete
	find . -type d -name '__pycache__' -delete

clean-test: ## remove test and coverage artifacts
	rm -rf .tox/ htmlcov/ .coverage

lint: ## check style with flake8
	flake8 geckopy tests

test: ## run tests quickly with the default Python
	pytest

test-all: ## run tests on every Python version with tox
	tox

coverage: ## check code coverage quickly with the default Python
		coverage run --source geckopy `which py.test`
		coverage report -m
		coverage html
		$(BROWSER) htmlcov/index.html

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/geckopy.rst docs/modules.rst
	sphinx-apidoc -o docs/ geckopy
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: clean copy-models ## package and upload a release
	python setup.py sdist bdist_wheel
	twine upload --skip-existing dist/*

dist: clean copy-models ## builds source and wheel package
	python setup.py sdist bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python setup.py install

copy-models: ## copy the latest model from main repository to package directory
	cp ../models/ecYeastGEM/ecYeastGEM.xml geckopy/data_files/ecYeastGEM_multi-pool.xml
	cp ../models/ecYeastGEM/ecYeastGEM_batch.xml geckopy/data_files/ecYeastGEM_single-pool.xml

protein-properties: ## generate the protein properties file from new downloads
	python fixtures/protein-properties.py
