# --------------------------------------------------
# Environments and top-level directories (no edits)
# --------------------------------------------------
# Canonical project root (works from any subdir)
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/.)

# Top-level directories
LOGS  := $(PROJECT_ROOT)/.logs
DATA  := $(PROJECT_ROOT)/data
R     := $(PROJECT_ROOT)/R
UTILS := $(R)/utils
OUT   := $(PROJECT_ROOT)/out
DRAFT := $(PROJECT_ROOT)/draft
FIGS  := $(PROJECT_ROOT)/figs

# Data
MAP_DATA   := $(OUT)/map-data.RData
NLOPT_DATA := $(OUT)/post-processed.RData

# --------------------------------------------------
# Check dependencies
# --------------------------------------------------
ifeq ($(shell command -v conda 2>/dev/null),)
  $(error !! ERROR: 'conda' not found in PATH)
endif

ifeq ($(shell conda info --envs | grep -q "$(CONDA_ENV)" && echo "exists"),)
  $(warning !! WARNING: Conda environment '$(CONDA_ENV)' not found)
  $(warning !!          Run 'make environments' to set up the workspace)
endif

ifeq ($(shell command -v Rscript 2>/dev/null),)
  $(error !! ERROR: 'Rscript' not found in PATH. This project requires R)
endif

# ifeq ($(shell command -v pandoc 2>/dev/null),)
#   $(error !! ERROR: 'pandoc' not found in PATH. Required for rendering the manuscript)
# endif

# --------------------------------------------------
# Logging
# --------------------------------------------------
DATE            := $(shell date +"%d-%m-%Y")
LOG_FILE        := $(LOGS)/log-$(DATE).log
LOGGER          := 2>&1 | tee -a $(LOG_FILE)
SUPPRESS_STDERR := 2>/dev/null
SUPPRESS_STDOUT := > /dev/null

# --------------------------------------------------
# Safe removal macro
# --------------------------------------------------
define SAFE_RM
	if [ -e "$(1)" ]; then \
	  ABS_PATH=$$(realpath "$(1)"); \
	  if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
	    echo " xx $$ABS_PATH"; \
	    rm -rf "$$ABS_PATH"; \
	  fi; \
	fi
endef

# --------------------------------------------------
# Loop removal macro
# --------------------------------------------------
define SAFE_RM_LIST
	for item in $(1); do $(call SAFE_RM,$$item); done
endef
