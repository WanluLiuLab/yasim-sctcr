.PHONY: dist
dist:
	python -m build

.PHONY: doc
doc:
	$(MAKE) -C doc

.PHONY: html
html:
	$(MAKE) -C doc html

.PHONY: pdf
pdf:
	$(MAKE) -C doc pdf

.PHONY: cleandoc
cleandoc:
	$(MAKE) -C doc clean
	$(MAKE) doc

.PHONY: serve-doc
serve-doc:
	python -m http.server -d doc/_build/html

.PHONY: pytype
pytype:
	 pytype --config=pytype.cfg src/labw_utils

.PHONY: test
test:
	$(MAKE) -c test
