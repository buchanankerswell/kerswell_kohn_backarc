# --------------------------------------------------
# Pandoc pipeline
# --------------------------------------------------
TEMPLATE      := $(DRAFT)/eisvogel.tex
TEMPLATE_NAME := $(notdir $(TEMPLATE))
BIB           := $(DRAFT)/main.bib

ifeq ($(TEMPLATE_NAME),agu.tex)
    CFG_YAML    := $(DRAFT)/agu.yaml
    CFG_SI_YAML := $(DRAFT)/agu-si.yaml
    CFG_SI_TPL  := $(DRAFT)/agu-si.tex
else ifeq ($(TEMPLATE_NAME),eisvogel.tex)
    CFG_YAML    := $(DRAFT)/eisvogel.yaml
    CFG_SI_YAML := $(DRAFT)/eisvogel-si.yaml
    CFG_SI_TPL  := $(TEMPLATE)
endif

BASE_PANDOC_FLAGS := --from markdown --filter pandoc-crossref $(PANDOC_FLAGS)

$(DRAFT)/manuscript.%: TPL                      := $(TEMPLATE)
$(DRAFT)/manuscript.%: META                     := $(CFG_YAML)
$(DRAFT)/manuscript.%: BIB_FLAG                 := --bibliography=$(BIB)
$(DRAFT)/manuscript.%: CFG_EXTRA                := --citeproc --csl=$(DRAFT)/agu.csl --pdf-engine=pdflatex

$(DRAFT)/supplementary-information.%: TPL       := $(CFG_SI_TPL)
$(DRAFT)/supplementary-information.%: META      := $(CFG_SI_YAML)
$(DRAFT)/supplementary-information.%: BIB_FLAG  := --bibliography=$(BIB)
$(DRAFT)/supplementary-information.%: CFG_EXTRA := --citeproc --csl=$(DRAFT)/agu.csl --pdf-engine=pdflatex

$(DRAFT)/diff-manuscript.%: META                := $(CFG_YAML)
$(DRAFT)/diff-supplementary-information.%: META := $(CFG_SI_YAML)

$(DRAFT)/response.%: TPL                        :=
$(DRAFT)/response.%: META                       := $(DRAFT)/response.yaml
$(DRAFT)/response.%: BIB_FLAG                   :=

TPL_FLAG     = $(if $(filter %.docx,$@),,$(if $(TPL),--template=$(TPL)))
META_FLAG    = $(if $(META),--metadata-file=$(META))
PANDOC_FLAGS = --standalone --from=markdown --filter=pandoc-crossref $(CFG_EXTRA) $(TPL_FLAG) $(BIB_FLAG) $(META_FLAG)

# Rendered files
DOCS := manuscript supplementary-information
PDFS := $(addprefix $(DRAFT)/, $(addsuffix .pdf, $(DOCS)))
DOCXS := $(addprefix $(DRAFT)/, $(addsuffix .docx, $(DOCS)))

# --------------------------------------------------
# Latexdiff pipeline
# --------------------------------------------------
# Configure pdflatex
PDFLATEX = pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT)

# Configure latexdiff
BACKUP_DATE ?= 15-jun-2026
DIFF_FLAGS   = --flatten --math-markup=0 --type=CFONT --graphics-markup=none --exclude-textcmd="ref,eqref,label" --add-to-config PICTUREENV=longtable --config="ARRENV=thebibliography|APACrefauthors"
