# Include definitions
include $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/../project.mk)
include local.mk

# Cleanup
EXTS       := aux bbl blg log toc out tex pdf
CLEAN      := $(foreach doc,$(DOCS) $(addprefix diff-,$(DOCS)),$(addprefix $(DRAFT)/$(doc).,$(EXTS))) $(addprefix $(DRAFT)/response.,$(EXTS))
DEEP_CLEAN :=

# Prevent make from deleting .tex files on failure
.PRECIOUS: $(DRAFT)/%.tex

# Targets
.PHONY: all manuscript diff clean deep-clean help

all: manuscript
manuscript: $(PDFS) $(DOCXS)
diff: all $(DRAFT)/diff-manuscript.pdf $(DRAFT)/diff-supplementary-information.pdf $(DRAFT)/response.pdf

$(DRAFT)/%.docx: $(DRAFT)/%.md
	@echo " -> $@"
	@$(if $(META),sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(META))
	@pandoc $< -o $@ $(PANDOC_FLAGS)

$(DRAFT)/%.pdf: $(DRAFT)/%.tex
	@echo " -> $@"
	@pandoc $(subst .tex,.md,$<) -o $@ $(PANDOC_FLAGS)

$(DRAFT)/bak/%.tex: ;

$(DRAFT)/%.tex: $(DRAFT)/%.md
	@echo " -> $@"
	@$(if $(META),sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(META))
	@pandoc $< -o $@ $(PANDOC_FLAGS)

$(DRAFT)/diff-%.pdf: $(DRAFT)/bak/%-$(BACKUP_DATE).tex $(DRAFT)/%.tex
	@echo " -> $@ [bak: $(BACKUP_DATE)]"
	@if [ ! -f "$<" ]; then \
	  echo "Error: Backup file $< not found!"; \
	  exit 1; \
	fi
	@latexdiff $(DIFF_FLAGS) $< $(word 2,$^) > $(DRAFT)/diff-$*.tex
	@$(PDFLATEX) $(DRAFT)/diff-$*.tex > /dev/null 2>&1
	@$(PDFLATEX) $(DRAFT)/diff-$*.tex > /dev/null 2>&1

clean:
	@$(call SAFE_RM_LIST,$(CLEAN))

deep-clean:
	@$(call SAFE_RM_LIST,$(DEEP_CLEAN))

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    manuscript Render manuscript and SI"
	@echo "    diff       Generate manuscript diff"
	@echo "    clean      Cleanup unnecessary files"
	@echo "    deep-clean Remove results and data (use with caution!)"
