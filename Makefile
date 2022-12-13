R = R/packages.R R/preprocess.R R/krige.R R/base-plots.R R/interpolation-plots.R R/summary-plots.R R/goutorbe-analysis.R R/download-data.R
DATAPURGE = data/hf.RData data/sectors.RData data/log*
DATACLEAN = data draft/assets/r
FIGSPURGE = figs draft/assets/figs

all: $(R)
	@./run.sh

purge:
	@rm -rf $(DATAPURGE) $(FIGSPURGE)

clean: purge
	@rm -rf $(DATACLEAN)

.PHONY: all purge clean purge