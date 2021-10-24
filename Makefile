R = packages.R preprocess.R krige.R base_plots.R interpolation_plots.R summary_plots.R
DATA = data/sa2006 data/tglobe
FIGS = figs/**/*.png
OPT = data/opt*.RData data/nloptr*

all: data/opt20.RData $(R) $(DATA)

data/opt20.RData:
	./run.sh

purge:
	@rm -f $(FIGS)

superclean: purge
	@rm -f $(OPT)

.PHONY: all clean