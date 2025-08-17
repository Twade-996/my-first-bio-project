# My First Bio Project
.hic

java -jar juicer_tools.jar dump observed KR ~/k562_test/rawdata/GSE63525_K562_combined_30.hic chr1 chr1 BP 100000 chr1_100kb.KRnorm.txt

chr1_100kb.KRnorm.txt

convert_3col_to_evr.py

chr1_100kb_evr_matrix.txt

evr.py

Chr1_100kb_evr_structure.pdb

----------------------------------------------------------------------------

.hic

java -Xmx16G -jar ~/software/juicer/juicer_tools.jar dump observed NONE ~/k562_test/rawdata/GSE63525_K562_combined_30.hic chr1 chr1 BP 100000 chr1_100kb_raw.txt

format_raw_hic.py

chr1_100kb.summary.formatted.txt & chr1_100kb_bins.bed

CscoreTool1.1 < windows.bed> < input.summary> < outputPrefix> < session> < minDis> [chrName]

Cscore(txt/bedgraph/…)

---------------------------------------------------------------------------

sample_1.bedgraph & sample_2.bedgraph

Compare_cscores.py

boundary_drift.png
compartment_switch_rate.png
cscore_compare.txt
cscore_comparison_summary.csv
cscore_delta_bar.png
cscore_distribution.png
cscore_heatmap.png
cscore_scatter.png
pearson.txt
spearman.txt

auto_plot_images.py

.png

----------------------------------------------------------------------------

sample_1.pdb & .bedgraph

python plot.py chr1_100kb_evr_structure.txt -cf chr1_100kb__cscore.txt --skip-range 1215 1425 --show

----------------------------------------------------------------------------

sample_1.pdb/.bedgraph & sample_2.pdb/.bedgraph

plot_interactive_morph.py & IO.py

----------------------------------------------------------------------------

sample_1.pdb/.bedgraph & sample_2.pdb/.bedgraph

plot_animation.py & IO.py

.mp4

