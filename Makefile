SHELL := /bin/bash
RSCRIPT ?= Rscript
R_BIN ?= R

.PHONY: help fast test check

help:
	@echo "Targets:"
	@echo "  fast   Run cheap validation (package tests only)"
	@echo "  test   Run unit tests via testthat"
	@echo "  check  Run R CMD check (no manual)"

# FAST should stay cheap and reliable; only runs unit tests.
fast: test

test:
	@$(RSCRIPT) -e 'if (!requireNamespace("testthat", quietly = TRUE)) { message("Install testthat to run tests"); quit(status = 2) }; testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)'

check:
	@$(R_BIN) CMD check --no-manual .
