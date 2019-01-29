# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from nearest_neighbor_threshold_analysis import generate_threshold_plot
from parse_sequences import txttoseqlist
import sys

# given a list of lists as input_repertoires (instead of files like above) run the threshold anaylsis
def run_threshold_analysis_repertoire_lists(output_directory, title, input_repertoires, k=1):
    all_seqs = []

    # read in the clusters
    repertoire_dict = {}
    i = 0
    for repertoire in input_repertoires:
        i += 1
        content = repertoire
        all_seqs += content
        for cdr in content:
            if cdr in repertoire_dict:
                found = False
                for nums in repertoire_dict[cdr]:
                    if nums == i:
                        found = True
                if not found:
                    repertoire_dict[cdr].append(i)
            else:
                repertoire_dict[cdr] = [i]

    num_clusters = i

    print(generate_threshold_plot(all_seqs, repertoire_dict, k, num_clusters, output_directory, title, step=.05))


# take in the 6 twin files
# combine each of the pairs of twins and the other 4 individuals into two separate repertoires
# call run_threshold_analysis 3 times, with each of the twin pairs separately
def twins_threshold_analysis(output_directory, twinA1, twinA2, twinC1, twinC2, twinD1, twinD2, beta=True):

    twin_pair_A = txttoseqlist(twinA1, trim_first=True) + txttoseqlist(twinA2, trim_first=True)
    twin_pair_C = txttoseqlist(twinC1, trim_first=True) + txttoseqlist(twinC2, trim_first=True)
    twin_pair_D = txttoseqlist(twinD1, trim_first=True) + txttoseqlist(twinD2, trim_first=True)

    # run threshold analysis on each twin pair
    if beta:
        run_threshold_analysis_repertoire_lists(output_directory, "TwA_Beta", [twin_pair_A, twin_pair_C + twin_pair_D], k=1)
        run_threshold_analysis_repertoire_lists(output_directory, "TwC_Beta", [twin_pair_C, twin_pair_A + twin_pair_D],  k=1)
        run_threshold_analysis_repertoire_lists(output_directory, "TwD_Beta", [twin_pair_D, twin_pair_A + twin_pair_C], k=1)
    else:
        run_threshold_analysis_repertoire_lists(output_directory, "TwA_Alpha", [twin_pair_A, twin_pair_C + twin_pair_D], k=1)
        run_threshold_analysis_repertoire_lists(output_directory, "TwC_Alpha", [twin_pair_C, twin_pair_A + twin_pair_D],  k=1)
        run_threshold_analysis_repertoire_lists(output_directory, "TwD_Alpha", [twin_pair_D, twin_pair_A + twin_pair_C], k=1)




# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 2:
    print("need the following arguments: output directory, repertoire files")
else:
    output_directory = sys.argv[1]
    repertoire_file_list = sys.argv[2].split(',')

twins_threshold_analysis(output_directory, repertoire_file_list[0], repertoire_file_list[1], repertoire_file_list[2],
                         repertoire_file_list[3], repertoire_file_list[4], repertoire_file_list[5], beta=False)

