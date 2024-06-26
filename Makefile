# Logging config
DATE = $(shell date +"%d-%m-%Y")
LOGFILE := log/log-$(DATE)
LOG := 2>&1 | tee -a $(LOGFILE)
# R scripts
R = R/nloptr.R \
		R/check-packages.R \
		R/download-assets.R \
		R/goutorbe-analysis.R \
		R/preprocess-map-data.R
# Directories with data and scripts
DATADIR = assets
# Cleanup directories
DATAPURGE = log
DATACLEAN = $(DATADIR)
FIGSPURGE = figs

all: nlopt

test: preprocess
	@R/test.R $(LOG)

nlopt: preprocess
	@R/nloptr.R $(LOG)

preprocess: $(DATADIR)
	@R/preprocess-map-data.R $(LOG)

$(DATADIR): check_packages
	@R/download-assets.R $(LOG);

check_packages: $(LOGFILE) $(R)
	@R/check-packages.R $(LOG)

$(LOGFILE):
	@if [ ! -e "$(LOGFILE)" ]; then \
		mkdir -p log; \
		touch $(LOGFILE); \
	fi

purge:
	@rm -rf $(DATAPURGE) $(FIGSPURGE)

clean: purge
	@rm -rf $(DATACLEAN)

.PHONY: clean purge check_packages preprocess nlopt test all