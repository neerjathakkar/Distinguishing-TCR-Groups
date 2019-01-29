# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

# command line arguments:
# for beta:
# python2 structural_cluster_histograms.py results/stcr_results/ STCR_global_beta_SW "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
# for alpha:
# python2 structural_cluster_histograms.py results/stcr_results/ STCR_global_alpha_SW "datasets/stcr_data/A3-9-A,datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A"


from matplotlib import pyplot as plt
import random


from sw_scoring import getDistanceSW
from parse_sequences import txttoseqlist
import seaborn as sns
import sys


# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 3:
    print("need the following arguments: output directory, file title, cluster files")
    sys.exit()
else:
    output_directory = sys.argv[1]
    title = sys.argv[2]
    cluster_file_list = sys.argv[3].split(',')


# return a list with the min scores of cdr1 compared to cdrs2, cdr_list1 = cdr_list2
def calculate_min_distances_within(cdr_list1, cdr_list2):
    scores = []
    # for each cdr in cdrs1, find the min score
    for cdr1 in cdr_list1:
        min = 1000000
        for cdr2 in cdr_list2:
            if cdr1 != cdr2:
                score = getDistanceSW(cdr1, cdr2, length_dep=True, gap_penalty=-10)
                if score < min:
                    min = score
        scores.append(min)
    return scores

# return a list with the min scores of cdr1 compared to cdrs2, cdr_list1 != cdr_list2
# same method as above except that the cdr in cdr list 2 must not be in cdr list 1
def calculate_min_distances_between(cdr_list1, cdr_list2):
    scores = []
    # for each cdr in cdrs1, find the min score
    for cdr1 in cdr_list1:
        min = 100000
        for cdr2 in cdr_list2:
            if cdr1 != cdr2 and cdr2 not in cdr_list1:
                score = getDistanceSW(cdr1, cdr2, length_dep=True, gap_penalty=-10)
                if score < min:
                    min = score
        scores.append(min)
    return scores


# generate a PDF and a CDF of the data
def generate_plots(self_scores, non_self_scores, output_directory, title):
    sns.set(font_scale=2)  # crazy big

    # density plot - PDF
    figure = sns.kdeplot(self_scores, shade=True)
    figure = sns.kdeplot(non_self_scores, shade=True)
    figure.set(ylim=(0, 5))
    fig = figure.get_figure()
    fig.show()
    fig.savefig(output_directory + title + "_kde_updated.pdf")

    plt.figure()    # new figure
    figure = sns.kdeplot(self_scores, shade=True, cumulative = True)
    figure = sns.kdeplot(non_self_scores, shade=True, cumulative = True)
    figure.set(ylim=(0, 1))
    fig = figure.get_figure()
    fig.savefig(output_directory + title + "_kde_cumulative_updated.pdf")

print("calculating self-scores...")
# calculate self-scores
all_self_scores = []
num_self_scores = 0
for file in cluster_file_list:
    cdrs = txttoseqlist(file)
    cdrs = set(cdrs)
    cdrs = list(cdrs)
    scores = calculate_min_distances_within(cdrs, cdrs)
    all_self_scores += scores
    num_self_scores += len(scores)
print(all_self_scores)

print("calculating non-self scores...")
# calculate non-self scores
all_non_self_scores = []
num_non_self_scores = 0
for file1 in cluster_file_list:
    for file2 in cluster_file_list:
        if file1 != file2:
            cdrs1 = list(set(txttoseqlist(file1)))
            cdrs2 = list(set(txttoseqlist(file2)))
            scores = calculate_min_distances_between(cdrs1, cdrs2)
            all_non_self_scores += scores
            num_non_self_scores += len(scores)

for i in range(len(all_non_self_scores)):
    if all_non_self_scores[i] == 1.0:
        all_non_self_scores[i] = random.uniform(0.9, 1.0)

print(all_non_self_scores)


generate_plots(all_self_scores, all_non_self_scores, output_directory, title)

