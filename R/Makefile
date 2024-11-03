# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
LOGS := $(PROJECT_ROOT)/logs
R := $(PROJECT_ROOT)/R
DATA := $(PROJECT_ROOT)/data
UTILS := $(R)/utils
OUT := $(PROJECT_ROOT)/out
FIGS := $(PROJECT_ROOT)/figs

# Logging
DATE := $(shell date +"%d-%m-%Y")
LOG_FILE := $(LOGS)/log-$(DATE).log
LOGGER := 2>&1 | tee -a $(LOG_FILE)

# Data
MAP_DATA := $(OUT)/map-data.RData
NLOPT_DATA := $(OUT)/nlopt-data.RData

# Cleanup targets
CLEAN := $(LOGS)/log*.log
DEEP_CLEAN := \
	      $(OUT)/map-data.RData \
	      $(OUT)/nlopt \
	      $(OUT)/nlopt-data.RData \
	      $(OUT)/heatflow \
	      $(OUT)/submap \
	      $(OUT)/relief \
	      $(FIGS)/global-dataset-composition.png \
	      $(FIGS)/sim-krg-accuracies.png \
	      $(FIGS)/nlopt.png \
	      $(FIGS)/submap

# Targets
.PHONY: all visualize test krige preprocess clean deep-clean help

all: preprocess krige visualize

visualize: $(LOG_FILE) $(MAP_DATA) $(NLOPT_data)
	@Rscript visualize.R $(UTILS) $(DATA) $(OUT) $(FIGS) $(LOGGER)

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
	@for item in $(CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
					echo " -> Safely removing: $$ABS_PATH"; \
					rm -rf "$$ABS_PATH"; \
				fi; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

deep-clean:
	@for item in $(DEEP_CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
					echo " -> Safely removing: $$ABS_PATH"; \
					rm -rf "$$ABS_PATH"; \
				fi; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    visualize         Visualize all results"
	@echo "    postprocess       Process kriging results"
	@echo "    krige             Krige with non-linear optimization"
	@echo "    preprocess        Fetch required datasets and preprocess"
	@echo "    clean             Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean        Deep clean figures and results (use with caution!)"
	@echo "    help              Show this help message"
