SHELL := /bin/bash

.PHONY: help fast test check

help:
	@echo "Targets:"
	@echo "  fast   Run cheap validation (package checks/tests)"
	@echo "  test   Run unit tests"
	@echo "  check  Run R CMD check"

# FAST should be cheap and reliable.
fast: test

test:
	@Rscript -e 'if (!requireNamespace("devtools", quietly=TRUE)) quit(status=2); devtools::test()'

check:
	@R CMD check .
