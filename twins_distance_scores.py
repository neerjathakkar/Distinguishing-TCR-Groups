# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from parse_sequences import txttoseqlist
from cluster_CDRs import cluster_data
from cluster_distance_scores import score_repertoire_twins, create_motif
import sys


def cluster_and_get_distance_scores(subject_file, subject_twin_file, other_files, clustering_distance, output_dir):
    log_file = open(output_dir + "output.txt", mode = 'w')
    log_file.write(subject_file[11:])

    Xlabels = txttoseqlist(subject_file, trim_first=True, remove_special=True)
    Xlabelstwin = txttoseqlist(subject_twin_file, trim_first=True, remove_special=True)

    print("clustering subject and their twin...")

    clusters = cluster_data(Xlabels, max_distance=clustering_distance)
    clusters_twin = cluster_data(Xlabelstwin, max_distance=clustering_distance)

    print("clustered subject and their twin")

    # convert clusters to list of clusters
    clusters_list = []
    for key in clusters:
        clusters_list.append(clusters[key])

    twin_clusters_list = []
    for key in clusters_twin:
        twin_clusters_list.append(clusters_twin[key])

    print("clustering other individuals...")

    other_clusterings = []
    # get the other twins reperotires
    for file in other_files:
        Xlabelsother = txttoseqlist(file, trim_first=True, remove_special=True)
        clusters2other = cluster_data(Xlabelsother, max_distance=clustering_distance)

        other_clusters = []
        for key in clusters2other:
            other_clusters.append(clusters2other[key])

        other_clusterings.append(other_clusters)

    print("clustered others")
    print("getting scores...")

    # get the scores
    scores = score_repertoire_twins(clusters_list,  twin_clusters_list, other_clusterings, min_clust_size=4)
    print(scores)

    n = 1

    for score_data in scores:
        log_file.write("\n\ncluster #" + str(n))
        # print(score_data)
        # print("parsing data...\n")

        # parse score data
        cluster = score_data[1]
        twin_cluster = score_data[0][1]
        twin_cluster_score = score_data[0][0]
        log_file.write("\ncluster is :" + str(cluster))
        log_file.write("\nthe closest twin cluster is " + str(twin_cluster))
        log_file.write("\ntheir score is: " + str(twin_cluster_score))

        # create motifs
        create_motif(cluster, str(n) + "_subject", output_dir)
        create_motif(twin_cluster, str(n) + "_subject_twin", output_dir)

        # calculate average distance to other clusters
        total_dist = 0
        count = 0
        other_scores = []

        for i in range(len(other_files)):
            data = score_data[i+2]
            repertoire = other_files[i][22:]
            other_score = data[0]
            other_cluster = data[2]
            log_file.write("\nscore with " + str(repertoire) + " is " + str(other_score))
            log_file.write("\ncluster: " + str(other_cluster))

            create_motif(other_cluster, str(n) + "_" + str(repertoire), output_dir)

            total_dist += other_score
            count += 1
            other_scores.append(other_score)

        avg_score = total_dist/count
        log_file.write("\naverage score to another cluster is " + str(avg_score))
        log_file.write("\nthe closest score to another cluster is " + str(min(other_scores)))

        log_file.write(' ')

        n += 1



# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 5:
    print("need the following arguments: output directory, threshold to cluster to, subject to analyze, their twin, all repertoire files")
    print("twin to analyze: TwA1, TwA2, TwC1, TwC2, TwD1, TwD2")
    sys.exit()

output_directory = sys.argv[1]
distance = float(sys.argv[2])
subject = sys.argv[3]
twin = sys.argv[4]
all_files = sys.argv[5].split(',')

filename = ""
twin_file = ""
for file in all_files:
    if subject in file:
        filename = file
    if twin in file:
        twin_file = file

other_files = []
for file in all_files:
    if subject not in file and twin not in file:
        other_files.append(file)

cluster_and_get_distance_scores(filename, twin_file, other_files, distance, output_directory)
