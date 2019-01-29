# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from sw_scoring import getDistanceSW

from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# takes two clusters C and C'
# returns the distance between them, using Smith Waterman
def between_cluster_distance(clust1, clust2):
    total_distance = 0
    count = 0

    for c1 in clust1:
        for c2 in clust2:
            total_distance += getDistanceSW(c1, c2, length_dep=True, gap_penalty=-10)
            count += 1

    return total_distance/float(count)


# takes a single cluster C and all of the clusters in another repertoire R
# returns the cluster C' in R with the closest distance to C, and that distance
def get_closest_cluster(cluster, other_repertoire, min_clust_size=2):
    min_distance = 1000
    closest_cluster = None

    for c in other_repertoire:
        if len(c) >= min_clust_size:
            dist = between_cluster_distance(cluster, c)
            if dist < min_distance:
                min_distance = dist
                closest_cluster = c

    return(closest_cluster, min_distance)


# repertoire: list of clusters
# other_repertoires: list of lists of clusters (the other repertoires, clustered)
# returns a sorted list with the closest cluster out of all the other repertoires
# use for Dash analysis
def score_repertoire(repertoire, other_repertoires):
    cluster_scores = []

    # score every cluster against each repertoire in the other_repertoires
    for cluster in repertoire:
        if len(cluster) >= 3:
            min_score = 1000
            partner_cluster = None
            closest_rep = 0

            rep_num = 1
            # find the closest cluster
            for other_rep in other_repertoires:
                (closest_cluster, score) = get_closest_cluster(cluster, other_rep)
                if score < min_score:
                    partner_cluster = closest_cluster
                    min_score = score
                    closest_rep = rep_num
                rep_num += 1

            # store the cluster, its pair cluster, which repertoire it came from, and their score
            cluster_scores.append((min_score, cluster, closest_rep, partner_cluster))

    # sort all of the results by score
    cluster_scores.sort()

    # return the sorted list
    return cluster_scores

# create a motif from a cluster
def create_motif(cluster, motif_name, out_dir):

    instances = []
    for seq in cluster:
        instances.append(Seq(seq, IUPAC.protein))

    m = motifs.create(instances, IUPAC.protein)
    file_name = out_dir + str(motif_name) + ".pdf"
    m.weblogo(file_name,
            show_xaxis=False,
            show_yaxis=False,
            show_errorbars=False,
            unit_name='',
            show_fineprint=False,
            format='pdf')


# return a list with the closest twin cluster and cluster in every other repertoire
def score_repertoire_twins(repertoire, twin_repertoire, other_repertoires, min_clust_size=3):
    cluster_scores = []
    i = 0

    # score every cluster against each repertoire in the other_repertoires
    for cluster in repertoire:
        if len(cluster) >= min_clust_size:

            # find the closest cluster in the twin
            (closest_cluster_twin, score_twin) = get_closest_cluster(cluster, twin_repertoire, min_clust_size=3)
            cluster_scores.append([(score_twin, closest_cluster_twin)])
            cluster_scores[i].append((cluster))

            j = 0
            # find the closest cluster in each of the other repertoires, store which repertoire we are in
            for other_rep in other_repertoires:
                (closest_cluster, score) = get_closest_cluster(cluster, other_rep, min_clust_size=3)
                # store the pair cluster, which repertoire it came from, and their score
                cluster_scores[i].append((score, j, closest_cluster))

                j += 1

            i += 1
    # sort all of the results by score
    cluster_scores.sort()

    # return the sorted list
    # contains the following: [(twin score, closest twin cluster), (the cluster itself),
    # (closest cluster and score for all other repertoires)]
    return cluster_scores
