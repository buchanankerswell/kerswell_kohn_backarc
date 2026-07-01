# Include definitions
include project.mk

# Targets
.PHONY: all build build-from-scratch build-from-repo manuscript visualize postprocess krige preprocess download upload environment clean deep-clean help

all: build

build: build-from-repo

build-from-scratch: environment preprocess krige postprocess visualize manuscript
	@open draft/manuscript.pdf

build-from-repo: environment download visualize manuscript
	@open draft/manuscript.pdf

manuscript:
	@$(MAKE) --no-print-directory -C $(DRAFT)

visualize: postprocess
	@$(MAKE) --no-print-directory -C $(R) visualize

postprocess: krige
	@$(MAKE) --no-print-directory -C $(R) postprocess

krige:
	@$(MAKE) --no-print-directory -C $(R) krige

preprocess:
	@$(MAKE) --no-print-directory -C $(R) preprocess

download:
	@$(MAKE) --no-print-directory -C $(R) download

upload: $(OUT)
	@$(MAKE) --no-print-directory -C $(R) upload

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
	@echo "    build-from-repo    Build study from data in OSF repo"
	@echo "    build-from-scratch Build study from scratch (time consuming)"
	@echo "    manuscript         Render manuscript"
	@echo "    visualize          Visualize results"
	@echo "    postprocess        Process Kriging results"
	@echo "    krige              Krige with non-linear optimization"
	@echo "    preprocess         Fetch and preprocess required datasets"
	@echo "    download           Download data from OSF repo"
	@echo "    environment        Create R environment"
	@echo "    clean              Remove generated files (safe)"
	@echo "    deep-clean         Remove results and data (use with caution!)"
