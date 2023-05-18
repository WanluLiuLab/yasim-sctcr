SPHINX_OPTS := -j auto

#SPHINX_APIDOC_ARGS := --separate \
#	--no-toc \
#	--module-first \
#	--maxdepth 2 \
#	--force \
#	--output-dir $(CURDIR)/_apidoc \
#	--doc-project 'API Docs'

BUILD_DIR := $(CURDIR)/_build
# SPHINX_APIDOC := PYTHONPATH="$(dir $(CURDIR))/src:${PYTHONPATH:-}" sphinx-apidoc

.PHONY: all
all: html

#.PHONY: apidoc
#apidoc:
#	rm -rf $(CURDIR)/_apidoc
#	$(SPHINX_APIDOC) $(SPHINX_APIDOC_ARGS) $(dir $(CURDIR))/src/labw_utils

.PHONY: preconfigure
preconfigure: src
	python preconfigure.py

.PHONY: src
src:
	$(MAKE) -C src

.PHONY: html
html: src preconfigure # apidoc  refs.bibtex.bib
	sphinx-build -M html $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

#refs.bibtex.bib: refs.bib
#	biber --tool \
#		--configfile=biberConf.xml \
#		--output-file refs.bibtex.bib \
#		--output-legacy-date \
#		--output-field-replace=location:address,journaltitle:journal \
#		refs.bib

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	$(MAKE) -C src clean

.PHONY: distclean
distclean: clean
