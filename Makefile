.PHONY: check, clean, flake8, wc, upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	pytest -v

clean:
	rm -fr _trial_temp .pytest_cache dist py3seq.egg-info
	find . -name '*.pyc' -print0 | $(XARGS) -0 rm
	find . -name '*~' -print0 | $(XARGS) -0 rm

flake8:
	find py3seq -name '*.py' -print0 | $(XARGS) -0 flake8

wc:
	find py3seq \( -name '*.py' -o -name '*.sh' \) -print0 | $(XARGS) -0 wc -l

# The upload target requires that you have access rights to PYPI. You'll
# also need twine installed (on OS X with brew, run 'brew install
# twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/py3seq-$$(grep __version__ py3seq/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
