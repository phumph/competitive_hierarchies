
manuscript_base=competition_v1
si_base=si

all: paper preprint si tables

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

clean:
	rm ${manuscript_base}.pdf ${manuscript_base}_ms.pdf ${si_base}.pdf ${si_base}_ms.pdf

.PHONY: all clean paper preprint tables growth-curves

# ANALYSIS SCRIPTS

growth-curves:
	Rscript analysis/scripts/fit_growth_traits.R \
		analysis/data/growthcurve_data_Pflu.txt \
		analysis/data/growthcurve_data_Psyr.txt \
		analysis/data/tsplits.csv \
		analysis/data/growth_traits_fitted.csv
		
