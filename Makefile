# Include definitions
include project.mk

# Targets
.PHONY: all build visualize postprocess test krige preprocess environment clean deep-clean help

all: build

build: environment preprocess krige postprocess visualize manuscript
	@open draft/manuscript.pdf

test:
	@$(MAKE) --no-print-directory -C $(R) test

download:
	@$(MAKE) --no-print-directory -C $(R) download

upload: $(OUT)
	@$(MAKE) --no-print-directory -C $(R) upload

visualize: postprocess
	@$(MAKE) --no-print-directory -C $(R) visualize

postprocess: krige
	@$(MAKE) --no-print-directory -C $(R) postprocess

krige:
	@$(MAKE) --no-print-directory -C $(R) krige

preprocess:
	@$(MAKE) --no-print-directory -C $(R) preprocess

manuscript:
	@$(MAKE) --no-print-directory -C $(DRAFT)

environment:
	@Rscript $(R)/environment.R

$(LOG_FILE):
	@if [ ! -e "$(LOG_FILE)" ]; then \
	  mkdir -p $(LOGS); \
	  touch $(LOG_FILE); \
	fi

clean:
	@$(MAKE) --no-print-directory -C $(R) clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) clean || true
	@find . -name ".DS_Store" -type f -delete

deep-clean: clean
	@$(MAKE) --no-print-directory -C $(R) deep-clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) deep-clean || true

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    visualize    Visualize results"
	@echo "    postprocess  Process Kriging results"
	@echo "    krige        Krige with non-linear optimization"
	@echo "    preprocess   Fetch and preprocess required datasets"
	@echo "    environment  Create R environment"
	@echo "    clean        Remove generated files (safe)"
	@echo "    deep-clean   Remove results and data (use with caution!)"
