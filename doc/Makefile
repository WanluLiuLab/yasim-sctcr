SPHINX_OPTS := -j auto

BUILD_DIR := $(CURDIR)/_build

.PHONY: all
all: html

.PHONY: preconfigure
preconfigure: src
	python preconfigure.py

.PHONY: src
src:
	$(MAKE) -C src

.PHONY: html
html: src preconfigure
	sphinx-build -M html $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) _cli_docs _root

.PHONY: distclean
distclean: clean
