R = packages.R preprocess.R krige.R base_plots.R interpolation_plots.R summary_plots.R
DATA = data/sa2006 data/tglobe data/hf.RData data/sectors.RData data/preprocess_log data/opt*.RData data/nloptr
FIGS = figs/*

all: data/opt20.RData $(R) $(DATA)

data/opt20.RData:
	@./run.sh

purge:
	@rm -rf $(FIGS)

clean: purge

superclean: purge
	@rm -rf $(DATA)

.PHONY: all clean