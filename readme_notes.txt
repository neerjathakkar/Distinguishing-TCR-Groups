Figure 2 (Oxford dataset):

PDFs/CDFs:
for beta:
python2 structural_cluster_histograms.py results/stcr_results/ STCR_global_beta "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
for alpha:
python2 structural_cluster_histograms.py results/stcr_results/ STCR_global_alpha "datasets/stcr_data/A3-9-A,datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A"

Nearest neighbor threshold plots:
for beta:
python2 run_threshold_analysis.py results/stcr_results/ STCR_global_beta "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
for alpha:
python2 run_threshold_analysis.py results/stcr_results/ STCR_global_alpha "datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A_only_10"

to generate closest within and between cluster pairs
python2 stcr_pairs.py results/stcr_results/ a3 "datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A_only_10"
python2 stcr_pairs.py results/stcr_results/ b3 "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"


Figure 3 (twins):

Nearest neighbor threshold plots:
python2 twins_threshold_analysis.py results/twin_results/ "datasets/twins_data/TwA1_B_top_1000.txt,datasets/twins_data/TwA2_B_top_1000.txt,datasets/twins_data/TwC1_B_top_1000.txt,datasets/twins_data/TwC2_B_top_1000.txt,datasets/twins_data/TwD1_B_top_1000.txt,datasets/twins_data/TwD2_B_top_1000.txt"
python2 twins_threshold_analysis.py results/twin_results/ "datasets/twins_data/TwA1_A_top_1000.txt,datasets/twins_data/TwA2_A_top_1000.txt,datasets/twins_data/TwC1_A_top_1000.txt,datasets/twins_data/TwC2_A_top_1000.txt,datasets/twins_data/TwD1_A_top_1000.txt,datasets/twins_data/TwD2_A_top_1000.txt"


Twins motifs:
python2 twins_distance_scores.py results/twin_results/motifs_distance_scoring/ 0.3 TwA1 TwA2
python2 twins_distance_scores.py results/twin_results/motifs_distance_scoring/dist_0.3_minsize_5_TwA1/ 0.3 TwA1 TwA2
python2 twins_distance_scores.py results/twin_results/motifs_distance_scoring/ 0.3 TwC1 TwC2

Figure 4 (Glanville):
HLA numbers (also generates background frequencies and analyzes their clusters)
python2 HLA_count_analysis.py results/glanville_results/ 0.2
python2 HLA_count_analysis.py results/glanville_results/ 0.3
python2 HLA_count_analysis.py results/glanville_results/ 0.4


Figure 5 (Dash):
Dash threshold plots:
python2 run_threshold_analysis.py results/dash_results/ Dash_beta_mice "datasets/dash_data/NP_beta.txt,datasets/dash_data/PA_beta.txt,datasets/dash_data/PB1_beta.txt,datasets/dash_data/PB1-F2_beta.txt"
python2 run_threshold_analysis.py results/dash_results/ Dash_beta_human "datasets/dash_data/BMLF1_beta.txt,datasets/dash_data/M1_beta.txt,datasets/dash_data/M38_beta.txt,datasets/dash_data/M45_beta.txt,datasets/dash_data/m139_beta.txt,datasets/dash_data/p65_beta.txt"

python2 run_threshold_analysis.py results/dash_results/ Dash_alpha_mice "datasets/dash_data/NP_alpha.txt,datasets/dash_data/PA_alpha.txt,datasets/dash_data/PB1_alpha.txt,datasets/dash_data/PB1-F2_alpha.txt"
python2 run_threshold_analysis.py results/dash_results/ Dash_alpha_human "datasets/dash_data/BMLF1_alpha.txt,datasets/dash_data/M1_alpha.txt,datasets/dash_data/M38_alpha.txt,datasets/dash_data/M45_alpha.txt,datasets/dash_data/m139_alpha.txt,datasets/dash_data/p65_alpha.txt"


Dash motifs:
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True BMLF1
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True p65
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True m139
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True M1

python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False M38
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False M45
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False NP
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PA
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PB1
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PB1-F2


Figure 6 (McPAS):
to get the results (and write to csv results file) for a pair of repertoires, call:
get_mcpas_results(FILE1, FILE2, output_csv_file)


