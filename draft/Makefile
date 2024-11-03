# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R := $(PROJECT_ROOT)/R
DRAFT := $(PROJECT_ROOT)/draft

# Manuscript
TEMPLATE := $(DRAFT)/eisvogel.tex
YAML := $(DRAFT)/eisvogel.yaml
SI_YAML := $(DRAFT)/eisvogel-si.yaml
# TEMPLATE := $(DRAFT)/agu.tex
# YAML := $(DRAFT)/agu.yaml
# SI_YAML := $(DRAFT)/agu-si.yaml
CSL := $(DRAFT)/agu.csl
BIB := $(DRAFT)/main.bib
DOC := $(DRAFT)/manuscript
SI := $(DRAFT)/supplementary-information
PREV_DOC := $(DRAFT)/bak/manuscript-28-oct-2025
DIFF_DOC := $(DRAFT)/diff

# Cleanup targets
CLEAN := $(DOC).{tex,pdf,aux,bbl,blg,log,toc} $(DIFF_DOC).{tex,pdf,aux,bbl,blg,log,toc} $(SI).{tex,pdf,aux,bbl,blg,log,toc}
# DEEP_CLEAN := $(PROJECT_ROOT)/figs
DEEP_CLEAN :=

# Targets
.PHONY: all diff manuscript supp update-date clean deep-clean help

all: manuscript

diff: $(DIFF_DOC).pdf

manuscript: $(DOC).tex $(DOC).pdf

supp: $(SI).tex $(SI).pdf

$(DIFF_DOC).pdf: $(PREV_DOC).tex $(DOC).tex
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript diff"
	@echo "    --------------------------------------------------"
	@echo " -> Comparing $(shell basename $(PREV_DOC).tex) with $(shell basename $(DOC).tex)"
	@latexdiff \
	  --flatten \
	  --math-markup=0 \
	  --type=CFONT \
	  --exclude-textcmd=caption \
	  --add-to-config PICTUREENV=longtable \
	  --graphics-markup=none \
	  $(PREV_DOC).tex $(DOC).tex > $(DIFF_DOC).tex 2> /dev/null
	@echo " -> Rendering $(shell basename $(DIFF_DOC).pdf) from $(shell basename $(DIFF_DOC).tex)"
	@pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(DIFF_DOC).tex > /dev/null 2>&1 || true
	@pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(DIFF_DOC).tex > /dev/null 2>&1 || true

$(DOC).pdf: update-date
	@echo " -> Rendering $(shell basename $(DOC).pdf) from $(shell basename $(DOC).tex)"
	@if [ "$(shell basename $(TEMPLATE))" = "eisvogel.tex" ]; then \
	  pandoc "$(DOC).md" -o "$(DOC).pdf" \
	    --from markdown \
	    --template "$(TEMPLATE)" \
	    --metadata-file=$(YAML) \
	    --bibliography=$(BIB) \
	    --csl=$(CSL) \
	    --filter=pandoc-crossref \
	    --citeproc \
	    --pdf-engine=pdflatex \
	    --pdf-engine-opt="-shell-escape" \
	    --pdf-engine-opt="-interaction=nonstopmode"; \
	else \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(DOC).tex > /dev/null 2>&1 || true; \
	  cd $(DRAFT) && bibtex $(notdir $(DOC)) > /dev/null 2>&1; \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(DOC).tex > /dev/null 2>&1 || true; \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(DOC).tex > /dev/null 2>&1 || true; \
	fi

$(SI).pdf: update-date
	@echo " -> Rendering $(shell basename $(SI).pdf) from $(shell basename $(SI).tex)"
	@if [ "$(shell basename $(TEMPLATE))" = "eisvogel.tex" ]; then \
	  pandoc "$(SI).md" -o "$(SI).pdf" \
	    --from markdown \
	    --template "$(TEMPLATE)" \
	    --metadata-file=$(SI_YAML) \
	    --bibliography=$(BIB) \
	    --csl=$(CSL) \
	    --filter=pandoc-crossref \
	    --citeproc \
	    --pdf-engine=pdflatex \
	    --pdf-engine-opt="-shell-escape" \
	    --pdf-engine-opt="-interaction=nonstopmode"; \
	else \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(SI).tex > /dev/null 2>&1 || true; \
	  cd $(DRAFT) && bibtex $(notdir $(SI)) > /dev/null 2>&1; \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(SI).tex > /dev/null 2>&1 || true; \
	  pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT) $(SI).tex > /dev/null 2>&1 || true; \
	fi

$(DOC).tex: update-date
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript"
	@echo "    --------------------------------------------------"
	@echo " -> Rendering $(shell basename $(DOC).md) with $(shell basename $(TEMPLATE)) template"
	@if [ "$(shell basename $(TEMPLATE))" = "eisvogel.tex" ]; then \
	  pandoc "$(DOC).md" -s -o "$(DOC).tex" \
	    --from markdown \
	    --template "$(TEMPLATE)" \
	    --metadata-file=$(YAML) \
	    --bibliography=$(BIB) \
	    --csl=$(CSL) \
	    --filter=pandoc-crossref \
	    --citeproc; \
	else \
	  pandoc "$(DOC).md" -o "$(DOC).tex" \
	    --from markdown \
	    --template "$(TEMPLATE)" \
	    --metadata-file=$(YAML) \
	    --bibliography=$(BIB) \
	    --filter=pandoc-crossref \
	    --natbib; \
	fi

$(SI).tex: update-date
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript"
	@echo "    --------------------------------------------------"
	@echo " -> Rendering $(shell basename $(SI).md) with $(shell basename $(TEMPLATE)) template"
	@if [ "$(shell basename $(TEMPLATE))" = "eisvogel.tex" ]; then \
	  pandoc "$(SI).md" -s -o "$(SI).tex" \
	    --from markdown \
	    --template "$(TEMPLATE)" \
	    --metadata-file=$(SI_YAML) \
	    --bibliography=$(BIB) \
	    --csl=$(CSL) \
	    --filter=pandoc-crossref \
	    --citeproc; \
	else \
	  pandoc "$(SI).md" -o "$(SI).tex" \
	    --from markdown \
	    --template agu-si.tex \
	    --metadata-file=$(SI_YAML) \
	    --bibliography=$(BIB) \
	    --csl=$(CSL) \
	    --filter=pandoc-crossref \
	    --natbib; \
	fi

update-date:
	@sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(DOC).md

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
	@echo "    clean       Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean  Deep clean figures and results (use with caution!)"
	@echo "    help        Show this help message"
