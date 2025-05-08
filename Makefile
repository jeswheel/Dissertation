# Variables
CHAPTER_DIRS = $(wildcard chapters/*)
RMW_FILES = $(foreach dir,$(CHAPTER_DIRS),$(wildcard $(dir)/*.Rnw))
TEX_FILES = $(RMW_FILES:.Rnw=.tex)
MAIN_TEX = main.tex
PDF_OUTPUT = main.pdf

# Default target
all: $(PDF_OUTPUT) dust

# Rule to create .tex files from .Rnw
%.tex: %.Rnw
	Rscript -e "knitr::knit('$<', output='$@')"

# Rule to build PDF from main.tex, including all dependencies
$(PDF_OUTPUT): $(MAIN_TEX) $(TEX_FILES)
	latexmk -pdf $(MAIN_TEX)

# Clean auxiliary files created by LaTeX
clean:
	latexmk -C
	rm -f $(TEX_FILES)

# Specifically clean only the output PDF
cleanpdf:
	rm -f $(PDF_OUTPUT)
	
# Remove only typical LaTeX auxiliary files
dust:
	rm -f $(wildcard *.log *.lof *.fls *.aux *.bbl *.blg *.fdb_latexmk *.lop *.out *.toc)

.PHONY: all clean cleanpdf dust