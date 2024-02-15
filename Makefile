# Logging config
DATE = $(shell date +"%d-%m-%Y")
LOGFILE := log/log-$(DATE)
LOG := 2>&1 | tee -a $(LOGFILE)
# R scripts
R = \
		R/check-packages.R \
		R/preprocess-map-data.R \
		R/preprocess-hf-data.R \
		R/krige.R \
		R/base-plots.R \
		R/interpolation-plots.R \
		R/summary-plots.R \
		R/goutorbe-analysis.R \
		R/download-assets.R
# Kriging parameters
MAXITR = 30
OPTALG = 2 # 1:Direct 2:SBPLX 3:NELDERMEAD 4:BOBYQA 5:COBYLA
INTERPWGHT = 0.5
VARIOGRAMWGHT = 0.5
NFOLD = 0.01 # 0:Leave-one-out 0<n<1: proportion of datapoints
NCORES = 6
# Directories with data and scripts
DATADIR = assets
# Cleanup directories
DATAPURGE = log assets.zip $(DATADIR)/hf.RData $(DATADIR)/sectors.RData
DATACLEAN = assets draft/assets/r/*.RData
FIGSPURGE = figs draft/assets/figs

all: preprocess
	@echo "=============================================" $(LOG)
	@./run.sh $(LOG)
	@echo "=============================================" $(LOG)

krige: preprocess
	@if [ ! -e "$(DATADIR)/hf_data/opt.RData" ]; then \
		R/krige $(MAXITR) $(OPTALG) $(INTERPWGHT) $(VARIOGRAMWGHT) $(NFOLD) $(NCORES) $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "Kriging interpolation found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

preprocess: $(LOGFILE) $(R) check_packages $(DATADIR)
	@if [ ! -e "$(DATADIR)/map_data/preprocessed-map-data.RData" ]; then \
		R/preprocess-map-data.R $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "Preprocessed map data found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi
	@if [ ! -e "$(DATADIR)/hf_data/preprocessed-hf-data.RData" ]; then \
		R/preprocess-hf-data.R $(LOG); \
		echo "=============================================" $(LOG); \
	else \
		echo "Preprocessed hf data found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

$(DATADIR): $(LOGFILE) $(R)
	@if [ ! -d "$(DATADIR)" ]; then \
		echo "=============================================" $(LOG); \
		R/download-assets.R $(LOG); \
	else \
		echo "Data assets found!" $(LOG); \
		echo "=============================================" $(LOG); \
	fi

check_packages: $(LOGFILE) $(R)
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