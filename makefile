
manuscript_base=competition_v1
si_base=si
submission="manuscript/competition_ismej.pdf"

all: paper preprint si tables figures

preprint:
	pandoc manuscript/meta_preprint.yaml \
		manuscript/${manuscript_base}.md \
		-o manuscript/${manuscript_base}.pdf \
		--template=manuscript/svm-latex-ms.tex \
		--filter pandoc-citeproc

si:
	pandoc manuscript/meta_si.yaml manuscript/${si_base}.md -o manuscript/${si_base}.pdf \
	-s \
	--template="manuscript/svm-latex-ms.tex" \
	--filter pandoc-citeproc

paper:
	pandoc manuscript/meta_paper.yaml \
		manuscript/${manuscript_base}.md \
		-o manuscript/${manuscript_base}_ms.pdf \
		-s \
		--template="manuscript/svm-latex-ms.tex" \
		--filter pandoc-citeproc

tables: manuscript/tables.tex
	bin/compile_tables analysis/tables manuscript/tables.tex

figures:
	pandoc manuscript/meta_figures.yaml \
		manuscript/figures.md \
		-o manuscript/figures.pdf \
		--template=manuscript/svm-latex-ms.tex \
		--filter pandoc-citeproc

submission:
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$(submission) manuscript/${manuscript_base}_ms.pdf manuscript/figures.pdf

clean:
	rm ${manuscript_base}.pdf ${manuscript_base}_ms.pdf ${si_base}.pdf ${si_base}_ms.pdf

.PHONY: all clean paper preprint tables submission figures si

