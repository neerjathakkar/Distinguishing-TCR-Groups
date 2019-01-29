# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information


# plots unidentified, correct, and incorrect CDR numbers using a basic KNN classifier
# only consider clusters > 3 unique CDRs considered

from sw_scoring import getDistanceSW
from parse_sequences import txttoseqlist

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

# calculate scores and return a dictionary of dictionaries of scores
def calculate_scores(cdr_list, length_dep=True):
    memo = {}
    scores = {}
    for cdr in cdr_list:
        cdr_dict = {}
        for other_cdr in cdr_list:
            if cdr != other_cdr:
                if (cdr, other_cdr) in memo:
                    cdr_dict[other_cdr] = memo[(cdr, other_cdr)]
                elif (other_cdr, cdr) in memo:
                    cdr_dict[other_cdr] = memo[(other_cdr, cdr)]
                else:
                    score = getDistanceSW(cdr, other_cdr, length_dep, gap_penalty=-10)
                    cdr_dict[other_cdr] = score
                    memo[(cdr, other_cdr)] = score
        scores[cdr] = cdr_dict
    return scores

# method that takes k, cdr, threshold, and the scores dict
# brute force returns the k closest neighbors
def k_closest_neighbors(k, cdr, threshold, scores_dict):
    found_neighbors = set()
    k_nearest = []

    for i in range(1, k+1):
        min = 2
        min_cdr = None
        for other_cdr in scores_dict:
            if not other_cdr in found_neighbors and not other_cdr == cdr:
                if scores_dict[cdr][other_cdr] < min:
                    min = scores_dict[cdr][other_cdr]
                    min_cdr = other_cdr
        if min <= threshold:
            k_nearest.append(min_cdr)
            found_neighbors.add(min_cdr)
        else:
            return None
    return k_nearest



def get_neighbors_clusters(k_nearest, structural_cluster_dict):
    clusters = []
    for cdr in k_nearest:
        if cdr in structural_cluster_dict:
            clusters.append(structural_cluster_dict[cdr])
    return clusters

# return majority element in a list
def majority(arr):
    c = Counter(arr)
    value, count = c.most_common()[0]
    return value


# takes a list of neighbor clusters - will be a list of lists of cluster numbers
# also takes the number of clusters
# ex. [0, 1], [2, 3], [1]
# if there is no majority cluster, return an empty list
def get_predicted_clusters(neighbor_clusters, num_clusters):
    # k value in KNN
    neighbors = len(neighbor_clusters)

    predicted = []
    twos = []
    for i in range(1, num_clusters+1):
        count = 0
        for arr in neighbor_clusters:
            if i in arr:
                count += 1
        if count == neighbors:
            predicted.append(i)
        if neighbors == 3 and count == 2:
            twos.append(i)
    if predicted == []:
        predicted = twos
    return predicted



# main function
# takes list of reads, k, step as parameters
# step is between 0 and 1 and says how to calculate t (ex. step = 0.2 -> t = 0, .2, .4, .6, .8, 1)
# produces a plot in dir
# cluster_number_to_rep: optional parameter that maps the cluster number to the repertoire it is from
def generate_threshold_plot(cdr_list, repertoire_dict, k, num_clusters, dir, file_title, step=.05,
         length_dep=True, confusion_matrix=False, cluster_number_to_rep=None):

    out_file = open(dir + file_title + "_log.txt", "w+")

    # calculate scores - all scores except self-scores, so do not need to worry about iterating over self
    scores = calculate_scores(cdr_list)

    # initialize a dictionary that maps a threshold value to the 3 counts
    threshold_totals = {}
    incorrect_cdrs = []
    incorrect_cdrs_set = set()
    incorrect_cluster_counter = {}

    # dict that keeps track of correctly IDed CDRs, their closest distance and CDRs, etc
    correctly_IDed_CDRs = {}

    if cluster_number_to_rep:
        out_file.write("how to map numbers to repertoires: ")
        out_file.write(str(cluster_number_to_rep))

    # iterate through thresholds
    for t in np.arange(0.0, 1.0 + step, step):
        print("step: " + str(t))

        # initialize lists/counts for unidentified, correct, and incorrect
        unidentified_count = 0
        correct_count = 0
        incorrect_count = 0


        # counts tuples of (actual cluster, predicted cluster)
        confusion_matrix_dict = {}


        # iterate through cdr_list
        for cdr in cdr_list:

            # find the k closest neighbors
            nearest_neighbors = k_closest_neighbors(k, cdr, t, scores)

             # if there are not k closest neighbors within the threshold -> unidentified
            if not nearest_neighbors:
                unidentified_count +=1
                continue

            neighbor_clusters = get_neighbors_clusters(nearest_neighbors, repertoire_dict)

            # get predicted clusters
            predicted_clusters = get_predicted_clusters(neighbor_clusters, num_clusters)

            if predicted_clusters == []:
                unidentified_count +=1
                continue

            matching_cluster = None

            # if this is correct -> correct
            incorrect_cluster = True
            for clust in predicted_clusters:
                if cdr in repertoire_dict:
                    for actual_clust in repertoire_dict[cdr]:
                        if clust == actual_clust:
                            incorrect_cluster = False
                            matching_cluster = clust

            if incorrect_cluster:
                for clust in predicted_clusters:
                    if cdr in repertoire_dict:
                        for actual_clust in repertoire_dict[cdr]:
                            if(actual_clust, clust) in confusion_matrix_dict:
                                confusion_matrix_dict[(actual_clust, clust)] += 1
                            else:
                                confusion_matrix_dict[(actual_clust, clust)] = 1

            else:
                if (matching_cluster, matching_cluster) in confusion_matrix_dict:
                    confusion_matrix_dict[(matching_cluster, matching_cluster)] += 1
                else:
                    confusion_matrix_dict[(matching_cluster, matching_cluster)] = 1

            if incorrect_cluster:
                incorrect_count += 1
                if cdr not in incorrect_cdrs_set:
                    if cdr in repertoire_dict:
                        incorrect_cdrs.append((cdr, repertoire_dict[cdr], nearest_neighbors, get_predicted_clusters(neighbor_clusters, num_clusters)))
                        incorrect_cdrs_set.add(cdr)
                        for clust in repertoire_dict[cdr]:
                            if clust in incorrect_cluster_counter:
                                incorrect_cluster_counter[clust] += 1
                            else:
                                incorrect_cluster_counter[clust] = 1
            else:
                correct_count += 1


        # populate threshold dictionary values
        threshold_totals[t] = (unidentified_count, correct_count, incorrect_count)

        print(t)
        print(confusion_matrix_dict)
        out_file.write(str(t) + "\n")
        out_file.write(str(confusion_matrix_dict))
        out_file.write("\n")

    # plot the threshold
    correct = []
    unidentified = []
    incorrect = []
    thresholds = []
    for key in sorted(threshold_totals.keys()):
        thresholds.append(key)
        correct.append(threshold_totals[key][1])
        incorrect.append(threshold_totals[key][2])
        unidentified.append(threshold_totals[key][0])


    plt.plot(thresholds, unidentified, '-o', color='g')
    plt.plot(thresholds, correct, '-o', color='b')
    plt.plot(thresholds, incorrect, '-o', color='m')

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)


    # title = file_title + " k: " + str(k)
    # plt.title(title)
    plt.savefig(dir + file_title + "_k_" + str(k) + ".pdf")
    print("saved fig")
    plt.figure()


    print("incorrectly identified")
    print(incorrect_cdrs)
    print(incorrect_cluster_counter)

    out_file.write("incorrectly identified\n")
    # out_file.write(str(incorrect_cdrs) + "\n")
    for incorrect_data in incorrect_cdrs:
        # print incorrect_data
        if cluster_number_to_rep:
            out_file.write(str(incorrect_data[0]) + " from " + str(cluster_number_to_rep[incorrect_data[1][0]][19:]) + \
                  " was incorrectly identified as being in  " + \
                  str(cluster_number_to_rep[incorrect_data[3][0]][19:]) + \
                " because its nearest neighbor is " + str(incorrect_data[2][0]) + "\n")
        else:
            out_file.write(str(incorrect_data[0]) + " from " + str(incorrect_data[1][0]) + \
                  " was incorrectly identified as being in  " + str(incorrect_data[3][0]) + \
                " because its nearest neighbor is " + str(incorrect_data[2][0]) + "\n")
    out_file.write(str(incorrect_cluster_counter) + "\n")

    out_file.write("totals: (unidentified, correct, incorrect) \n")
    out_file.write(str(threshold_totals))

    out_file.close()

    return threshold_totals




# given a list of files input_repertoires read them in and run the threshold anaylsis
def run_threshold_analysis(output_directory, title, input_repertoires, trim_first=True, rawtext=True, k=1):
    all_seqs = []

    # read in the clusters
    repertoire_dict = {}
    cluster_number_to_rep = {}
    i = 0
    for file in input_repertoires:
        i += 1
        cluster_number_to_rep[i] = file
        content = txttoseqlist(file, trim_first=trim_first, rawtext=rawtext)
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

        # get cluster length
        set_content = set(content)
    num_clusters = i

    print(generate_threshold_plot(all_seqs, repertoire_dict, k, num_clusters, output_directory, title, step=.05, cluster_number_to_rep=cluster_number_to_rep))


# returns (unidentified, correct, incorrect) at the list of specified thresholds
def threshold_analysis_no_plot(cdr_list, repertoire_dict, k, num_clusters, step=.05):
    # calculate scores - all scores except self-scores, so do not need to worry about iterating over self
    scores = calculate_scores(cdr_list)

    # initialize a dictionary that maps a threshold value to the 3 counts
    threshold_totals = {}

    # iterate through thresholds
    t = 0
    while t <= 1.0:
        # initialize lists/counts for unidentified, correct, and incorrect
        unidentified_count = 0
        correct_count = 0
        incorrect_count = 0


        # iterate through cdr_list
        for cdr in cdr_list:

            # find the k closest neighbors
            nearest_neighbors = k_closest_neighbors(k, cdr, t, scores)

             # if there are not k closest neighbors within the threshold -> unidentified
            if not nearest_neighbors:
                unidentified_count +=1
                continue

            neighbor_clusters = get_neighbors_clusters(nearest_neighbors, repertoire_dict)

            # get predicted clusters
            predicted_clusters = get_predicted_clusters(neighbor_clusters, num_clusters)

            if predicted_clusters == []:
                unidentified_count +=1
                continue

            # if this is correct -> correct
            incorrect_cluster = True
            for clust in predicted_clusters:
                if cdr in repertoire_dict:
                    for actual_clust in repertoire_dict[cdr]:
                        if clust == actual_clust:
                            incorrect_cluster = False


            if incorrect_cluster:
                incorrect_count += 1
            else:
                correct_count += 1


        # populate threshold dictionary values
        threshold_totals[t] = (unidentified_count, correct_count, incorrect_count)

        t = round(t + step, 2)

    return threshold_totals


# given a list of files input_repertoires read them in and run the threshold anaylsis
# use for McPAS results
def run_threshold_analysis_no_plot(input_repertoires, pathology_dict, k=1):
    all_seqs = []

    # read in the clusters
    repertoire_dict = {}
    i = 0
    for rep in input_repertoires:
        i += 1
        print(rep)
        content = pathology_dict[rep]
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

    return threshold_analysis_no_plot(all_seqs, repertoire_dict, k, num_clusters, step=.05)
