# # Latex Makefile using latexmk
# # Modified by Dogukan Cagatay <dcagatay@gmail.com>
# # Originally from : http://tex.stackexchange.com/a/40759
# # Modified for compilation of Markdown files with ideas from
# # http://plain-text.co/pull-it-together.html
#
# ## What extension (e.g. md, markdown, mdown) is being used
# ## for markdown files MEXT = md
# MEXT = md
#
# ## Expands to a list of all markdown files in the working directory
# SRC = $(wildcard *.$(MEXT))
#
# ## working bibliography file (make sure there is only one)
# BIB = $(wildcard *.bib)
#
# ## Target file names with extensions (use substitution references)
# ## x.pdf depends on x.md, x.html depends on x.md, etc
# ## https://www.gnu.org/software/make/manual/make.html#Substitution-Refs
# ## pdfs can be built from all md and tex files
# ## the sort function removes duplicates
# PDF_TEX = $(SRC:.tex=.pdf)
# PDF_MD = $(SRC:.md=.pdf)
# PDFS = $(sort $(PDF_MD) $(PDF_TEX))
# HTML = $(SRC:.md=.html)
# TEX = $(SRC:.md=.tex)
# DOCX = $(SRC:.md=.docx)
#
# ifeq ($(strip $(BIB)),)
# else
# 	PANDOC += --filter pandoc-citeproc --bibliography=$(BIB)
# endif
#
# ## Compile commands
#
# # MAIN LATEXMK RULE
# # -pdf tells latexmk to generate PDF directly (instead of DVI).
# # -pdflatex="" tells latexmk to call a specific backend with specific options.
# # -use-make tells latexmk to call make for generating missing files.
# # -interactive=nonstopmode keeps the pdflatex backend from stopping at a
# # missing file reference and interactively asking you for an alternative.
# LATEXMK = latexmk -use-make -pdf -pdflatex="pdflatex -interactive=nonstopmode"
#
# PANDOC  = pandoc -V papersize:a4 -V urlcolor:blue \
# 	--from markdown+autolink_bare_uris+simple_tables+table_captions+yaml_metadata_block
#
#
# 	# You want latexmk to *always* run, because make does not have all the info.
# 	# Also, include non-file targets in .PHONY so they are run regardless of any
# 	# file of the given name existing.
# 	.PHONY: $(PDF_TEX) all clean bib
#
#
# ## Rules
# all: clean bib $(PDFS)
# pdf:    clean bib $(PDFS)
# html:   $(HTML)
# tex:    $(TEX)
# docx:   $(DOCX)
#
# ## Implicit rules
# ## https://www.gnu.org/software/make/manual/html_node/Using-Implicit.html#Using-Implicit
# ## Each implicit rule has a target pattern and prerequisite patterns.
# ## There may be many implicit rules with the same target pattern
# ## The rule that actually applies is the one whose prerequisites exist or can be made
#
# ## '$@' is a variable holding the name of the target,
# ## and '$<' is a variable holding the (first) dependency of a rule.
#
# %.pdf: %.md
# 	$(PANDOC) -o $@ $<
#
# %.pdf: %.tex
# 	$(LATEXMK) $<
#
# cleanall:
# 	latexmk -C
#
# clean:
# 	-latexmk -c
#

manuscript_file=competition_v1.md
si_file=si.md
bib_file=competition.bib

all: competition_ms.md competition_ms.tex competition_ms.pdf

competition_ms.md: ${manuscript_file} ${bib_file}
		pandoc "${manuscript_file}" -o competition_ms.md \
	    -f html \
			-t markdown \
			-s \
			--mathjax \
			--standalone \
			--bibliography="${bib_file}" \
			--csl=evolution.csl \
			--css=pandoc.css \
			--filter=pandoc-citeproc

competition_ms.tex: competition_ms.md
	pandoc competition_ms.md -o competition_ms.tex \
		  -f html \
			-t latex \
			--pdf-engine=xelatex \
			--include-in-header=format.sty

competition_ms.pdf: competition_ms.tex
	xelatex competition_ms.tex competition_ms.pdf
