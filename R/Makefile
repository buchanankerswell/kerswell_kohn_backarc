# Include definitions
include $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/../project.mk)

# Cleanup
CLEAN := $(LOGS)/log*.log
DEEP_CLEAN := $(OUT) $(FIGS)

# Targets
.PHONY: all download upload visualize postprocess test krige preprocess clean deep-clean help

all: preprocess krige postprocess visualize

download: $(LOG_FILE) $(MAP_DATA) $(NLOPT_data)
	@Rscript osf.R $(UTILS) $(DATA) "download" $(LOGGER)
	@Rscript osf.R $(UTILS) $(OUT) "download" $(LOGGER)

upload: $(LOG_FILE) $(MAP_DATA) $(NLOPT_data)
	@Rscript osf.R $(UTILS) $(DATA) "upload" $(LOGGER)
	@Rscript osf.R $(UTILS) $(OUT) "upload" $(LOGGER)

visualize: $(LOG_FILE) $(MAP_DATA) $(NLOPT_data)
	@Rscript visualize.R $(UTILS) $(DATA) $(OUT) $(FIGS) $(LOGGER)

postprocess: $(LOG_FILE) $(MAP_DATA) $(NLOPT_data)
	@Rscript postprocess.R $(UTILS) $(DATA) $(OUT) $(FIGS) $(LOGGER)

test: $(LOG_FILE) $(MAP_DATA)
	@Rscript test.R $(UTILS) $(DATA) $(OUT) $(FIGS) $(LOGGER)

krige: $(LOG_FILE) $(MAP_DATA)
	@Rscript krige.R $(UTILS) $(DATA) $(OUT) $(FIGS) $(LOGGER)

preprocess: $(LOG_FILE)
	@Rscript preprocess.R $(UTILS) $(DATA) $(OUT) $(LOGGER)

$(LOG_FILE):
	@if [ ! -e "$(LOG_FILE)" ]; then \
	  mkdir -p $(LOGS); \
	  touch $(LOG_FILE); \
	fi

clean:
	@$(call SAFE_RM_LIST,$(CLEAN))

deep-clean:
	@$(call SAFE_RM_LIST,$(DEEP_CLEAN))

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    download    Download data from OSF repo"
	@echo "    upload      Upload data from OSF repo"
	@echo "    visualize   Visualize all results"
	@echo "    postprocess Process interpolation results"
	@echo "    krige       Krige with non-linear optimization"
	@echo "    preprocess  Fetch required datasets and preprocess"
	@echo "    clean       Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean  Deep clean figures and results (use with caution!)"
	@echo "    help        Show this help message"
