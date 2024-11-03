# Top-level dirs
PROJECT_ROOT := $(CURDIR)
R := $(PROJECT_ROOT)/R
DRAFT := $(PROJECT_ROOT)/draft

# Targets
.PHONY: build all visualize postprocess test krige preprocess environment clean deep-clean help

all: build

build: check-deps environment preprocess krige visualize manuscript
	@echo "    --------------------------------------------------"
	@echo "    Study built successfully!"
	@echo "    --------------------------------------------------"
	@open draft/manuscript.pdf

visualize:
	@$(MAKE) --no-print-directory -C $(R) visualize

postprocess: krige
	@$(MAKE) --no-print-directory -C $(R) postprocess

test:
	@$(MAKE) --no-print-directory -C $(R) test

krige:
	@$(MAKE) --no-print-directory -C $(R) krige

preprocess: $(LOG_FILE)
	@$(MAKE) --no-print-directory -C $(R) preprocess

manuscript:
	@$(MAKE) --no-print-directory -C $(DRAFT)

environment: $(LOG_FILE)
	@Rscript $(R)/environment.R

check-deps:
	@echo "    --------------------------------------------------"
	@echo "    Checking required dependencies"
	@echo "    --------------------------------------------------"
	@if ! command -v R >/dev/null 2>&1; then \
	    echo "ERROR: R not found in PATH."; \
	    exit 1; \
	fi
	@if ! command -v pandoc >/dev/null 2>&1; then \
	    echo " !! ERROR: pandoc not found in PATH!"; \
	    exit 1; \
	fi
	@echo " -- All dependencies found"

clean:
	@echo "    --------------------------------------------------"
	@echo "    Cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(R) clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) clean || true
	@find . -name ".DS_Store" -type f -delete

deep-clean: clean
	@echo "    --------------------------------------------------"
	@echo "    Deep cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(R) deep-clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) deep-clean || true

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    visualize         Visualize all results"
	@echo "    postprocess       Process kriging results"
	@echo "    krige             Krige with non-linear optimization"
	@echo "    preprocess        Fetch required datasets and preprocess"
	@echo "    environment       Check and install R packages"
	@echo "    clean             Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean        Deep clean results, figures, and data files (use with caution!)"
	@echo "    help              Show this help message"
