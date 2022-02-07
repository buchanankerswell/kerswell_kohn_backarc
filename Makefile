R = R/packages.R R/preprocess.R R/krige.R R/base_plots.R R/interpolation_plots.R R/summary_plots.R
DATA = data/sa2006 data/tglobe data/utig
DATAPURGE = data/hf.RData data/sectors.RData data/preprocess_log
DATACLEAN = data/opt*.RData data/nloptr
FIGS = figs/*

all: data/opt20.RData $(R) $(DATA)

data/opt20.RData:
	@./run.sh

purge:
	@rm -rf $(FIGS) $(DATAPURGE)

clean: purge

superclean: purge
	@rm -rf $(DATACLEAN)

.PHONY: all clean purge superclean data/opt20.RData