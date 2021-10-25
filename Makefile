R = packages.R preprocess.R krige.R base_plots.R interpolation_plots.R summary_plots.R
DATA = data/sa2006 data/tglobe
HF = data/hf.RData
FIGS = figs/*
LOGS = data/preprocess_log
OPT = data/opt*.RData data/nloptr

all: data/opt20.RData $(R) $(DATA)

data/opt20.RData:
	@./run.sh

purge:
	@rm -rf $(FIGS)

clean: purge

superclean: purge
	@rm -rf $(OPT) $(HF) $(LOGS)

.PHONY: all clean