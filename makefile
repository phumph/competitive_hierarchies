
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

.PHONY: all preprint si paper tables clean


# ANALYSIS SCRIPTS

SRC=analysis/scripts
DAT=analysis/data
FIG=analysis/figs
TAB=analysis/tables

.PHONY: growth-curves comp-traits pca pairwise traits interactions analysis facilitations analysis

analysis: traits interactions facilitations

traits: growth-traits comp-traits pca pairwise

growth-traits: $(DAT)/growth_traits_fitted.csv $(FIG)/growth_curves.pdf

# files that are NOT GENERATED:
# $(DAT)/growthcurve_data_Pflu.txt
# $(DAT)/growthcurve_data_Psyr.txt 
# $(DAT)/tsplits.csv
# $(DAT)/c_matrix.txt 
# $(DAT)/i_matrix.txt
# $(DAT)/strain_metadata.csv
# $(DAT)/phylogenetic_distance.txt

$(DAT)/growth_traits_fitted.csv: $(SRC)/fit_growth_traits.R $(DAT)/growthcurve_data_Pflu.txt $(DAT)/growthcurve_data_Psyr.txt $(DAT)/tsplits.csv
	Rscript --vanilla $< \
		$(DAT)/growthcurve_data_Pflu.txt \
		$(DAT)/growthcurve_data_Psyr.txt \
		$(DAT)/tsplits.csv \
		$@
			
$(FIG)/growth_curves.pdf: $(SRC)/make_growth_curves.R $(DAT)/growthcurve_data_Pflu.txt $(DAT)/growthcurve_data_Psyr.txt
	Rscript --vanilla $< \
		$(DAT)/growthcurve_data_Pflu.txt \
		$(DAT)/growthcurve_data_Psyr.txt \
		$@

comp-traits: $(DAT)/comp_traits.csv

$(DAT)/comp_traits.csv: $(SRC)/calc_comp_traits.R $(DAT)/c_matrix.txt $(DAT)/i_matrix.txt
	Rscript --vanilla $^ $@

pca: $(DAT)/all_traits.txt $(DAT)/pca_traits.txt $(TAB)/trait_summary_stats.tex

$(DAT)/all_traits.txt $(DAT)/pca_traits.txt $(TAB)/trait_summary_stats.tex: $(SRC)/make_traits_analysis.R $(DAT)/comp_traits.csv $(DAT)/growth_traits_fitted.csv $(DAT)/strain_metadata.csv
	Rscript --vanilla $^ $(DAT)/all_traits.txt $(DAT)/pca_traits.txt $(TAB)/trait_summary_stats.tex

pairwise: $(FIG)/pairwise_traits_biplots_comp.pdf $(FIG)/pairwise_traits_biplots_growth.pdf $(FIG)/pairwise_traits_biplots_growth_comp.pdf

$(FIG)/pairwise_traits_biplots_comp.pdf $(FIG)/pairwise_traits_biplots_growth.pdf $(FIG)/pairwise_traits_biplots_growth_comp.pdf: $(SRC)/make_pairwise_trait_plots.R $(DAT)/all_traits.txt
	Rscript --vanilla $^ $(FIG)/pairwise_traits_biplots

interactions: $(FIG)/mv_trait_dists.png $(FIG)/interaction_barplot.pdf $(TAB)/mv_dist_res.tex $(TAB)/lm_trait-v-phylo-dist_res.tex $(TAB)/mn_outcomes_res.tex $(TAB)/mn_coef_res.tex

$(DAT)/interaction_pairs.txt: $(SRC)/make_interaction_pairs.R $(DAT)/c_matrix.txt $(DAT)/strain_metadata.csv $(DAT)/phylogenetic_distance.txt
	Rscript --vanilla $^ $@

$(FIG)/mv_trait_dists.png $(FIG)/interaction_barplot.pdf $(TAB)/mv_dist_res.tex $(TAB)/lm_trait-v-phylo-dist_res.tex $(TAB)/mn_outcomes_res.tex $(TAB)/mn_coef_res.tex: $(SRC)/make_outcomes_analysis.R $(DAT)/all_traits.txt $(DAT)/interaction_pairs.txt $(DAT)/pca_traits.txt
	Rscript --vanilla $^

facilitations: $(FIG)/facil_plot_summaries.pdf $(FIG)/facil_plot_Psyr_ranks.pdf $(FIG)/facil_plot_Pflu_ranks.pdf $(DAT)/facil_effects_summary.csv

$(FIG)/facil_plot_summaries.pdf $(FIG)/facil_plot_Psyr_ranks.pdf $(FIG)/facil_plot_Pflu_ranks.pdf $(DAT)/facil_effects_summary.csv: $(SRC)/calc_facil_effects.R $(DAT)/c_matrix.txt $(DAT)/i_matrix.txt $(DAT)/all_traits.txt
	Rscript --vanilla $^ $(DAT)/facil_effects_summary.csv $(FIG)

clean-analysis:
	rm -f $(FIG)/facil_plot_summaries.pdf $(FIG)/facil_plot_Psyr_ranks.pdf $(FIG)/facil_plot_Pflu_ranks.pdf $(DAT)/facil_effects_summary.csv
	rm -f $(FIG)/mv_trait_dists.png $(FIG)/interaction_barplot.pdf $(TAB)/mv_dist_res.tex $(TAB)/lm_trait-v-phylo-dist_res.tex $(TAB)/mn_outcomes_res.tex $(TAB)/mn_coef_res.tex
	rm -f $(DAT)/interaction_pairs.txt
	rm -f $(FIG)/pairwise_traits_biplots_comp.pdf $(FIG)/pairwise_traits_biplots_growth.pdf $(FIG)/pairwise_traits_biplots_growth_comp.pdf
	rm -f $(DAT)/all_traits.txt $(DAT)/pca_traits.txt $(TAB)/trait_summary_stats.tex
	rm -f $(DAT)/comp_traits.csv
	rm -f $(FIG)/growth_curves.pdf
	rm -f $(DAT)/growth_traits_fitted.csv
