# Logging config
DATE = $(shell date +"%d-%m-%Y")
LOGFILE := log/log-$(DATE)
LOG := 2>&1 | tee -a $(LOGFILE)
# R scripts
R = \
		R/check-packages.R \
		R/preprocess-map-data.R \
		R/nloptr.R \
		R/interpolation-plots.R \
		R/summary-plots.R \
		R/goutorbe-analysis.R \
		R/download-assets.R
# Directories with data and scripts
DATADIR = assets
# Cleanup directories
DATAPURGE = log
DATACLEAN = $(DATADIR)
FIGSPURGE = figs

all: krige

nlopt: preprocess
	@if [ ! -e "$(DATADIR)/nlopt_data/interpolation-summary.RData" ]; then \
		R/nloptr.R $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "nlopt results found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

preprocess: $(DATADIR)
	@if [ ! -e "$(DATADIR)/map_data/map-data.RData" ]; then \
		R/preprocess-map-data.R $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "Preprocessed map data found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

$(DATADIR): check_packages
	@if [ ! -d "$(DATADIR)" ]; then \
		R/download-assets.R $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "Data assets found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

check_packages: $(LOGFILE) $(R)
	@echo "=============================================" $(LOG)
	@R/check-packages.R $(LOG)
	@echo "=============================================" $(LOG)

$(LOGFILE):
	@if [ ! -e "$(LOGFILE)" ]; then \
		mkdir -p log; \
		touch $(LOGFILE); \
	fi

purge:
	@rm -rf $(DATAPURGE) $(FIGSPURGE)

clean: purge
	@rm -rf $(DATACLEAN)

.PHONY: clean purge check_packages preprocess all