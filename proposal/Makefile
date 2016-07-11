OUTPUT_DIR=out
DOCUMENTS=$(wildcard *.tex)
PDFS=$(addprefix $(OUTPUT_DIR)/,$(DOCUMENTS:.tex=.pdf))

LATEXMK=latexmk
TEX_OPTIONS=-pdf -shell-escape

GRAPHICS_DIR=graphics
METAPOST=mpost
INKSCAPE=inkscape

MP_GRAPHICS=$(wildcard $(GRAPHICS_DIR)/*.mp)
MP_GRAPHICS_1=$(subst $(GRAPHICS_DIR),$(OUTPUT_DIR),$(MP_GRAPHICS:.mp=.1))
SVG_GRAPHICS=$(wildcard $(GRAPHICS_DIR)/*.svg)
SVG_GRAPHICS_EPS=$(subst $(GRAPHICS_DIR),$(OUTPUT_DIR),$(SVG_GRAPHICS:.svg=.eps))

all: out_dir mp_graphics $(PDFS)

out_dir:
	mkdir -p $(OUTPUT_DIR)

clean:
	- rm -rf $(OUTPUT_DIR)

$(OUTPUT_DIR)/%.pdf: %.tex
	$(LATEXMK) $(TEX_OPTIONS) -jobname=$(subst .pdf,,$@) $<

mp_graphics: out_dir $(MP_GRAPHICS_1)

$(OUTPUT_DIR)/%.1: $(GRAPHICS_DIR)/%.mp
	$(METAPOST) --jobname $(subst .1,,$@) $<
	mv $(subst .1,*,$(subst $(OUTPUT_DIR)/,,$@)) $(OUTPUT_DIR)

svg_graphics: out_dir $(SVG_GRAPHICS_EPS)

$(OUTPUT_DIR)/%.eps: $(GRAPHICS_DIR)/%.svg
	$(INKSCAPE) --file $< --export-eps=$@

watch: all
	when-changed *.tex *.bib -c make
