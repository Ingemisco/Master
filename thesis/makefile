TIKZ_SRC_DIR := tikz-tex
TIKZ_OUT_DIR := tikz-fig

MAIN_TEX     := main.tex 

TIKZ_TEX     := $(wildcard $(TIKZ_SRC_DIR)/*.tex)
TIKZ_PDF     := $(patsubst $(TIKZ_SRC_DIR)/%.tex,$(TIKZ_OUT_DIR)/%.pdf,$(TIKZ_TEX))

all: $(TIKZ_TEX) main

$(TIKZ_OUT_DIR)/%.pdf: $(TIKZ_SRC_DIR)/%.tex | $(TIKZ_OUT_DIR)
	pdflatex -output-directory=$(TIKZ_OUT_DIR) -interaction=nonstopmode $<

$(TIKZ_OUT_DIR):
	mkdir -p $(TIKZ_OUT_DIR)

main: $(TIKZ_PDF)
	latexmk -pdf $(MAIN_TEX)

clean:
	rm -r $(TIKZ_OUT_DIR)
	latexmk -c $(MAIN_TEX)


