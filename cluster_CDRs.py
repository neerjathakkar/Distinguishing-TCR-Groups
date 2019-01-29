# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqIO.FastaIO import SimpleFastaParser

from sw_scoring import getDistanceSW

def construct_distance_vector(Xlabels, length_dep=True, gap_penalty=-10):
    distances = []
    for i in range(len(Xlabels)):
        for j in range(i+1, len(Xlabels)):
            dist = getDistanceSW(Xlabels[i], Xlabels[j], length_dep=length_dep, gap_penalty=gap_penalty)
            distances.append(dist)
    return distances

def clusters_to_dict(cluster_arr, Xlabels):
    clusters = {}
    for i in range(len(cluster_arr)):
        if cluster_arr[i] in clusters:
            clusters[cluster_arr[i]].append(Xlabels[i])
        else:
            clusters[cluster_arr[i]] = [Xlabels[i]]
    return clusters


def cluster_data(Xlabels, max_distance=0.2, colors_dict=False, structural_cluster_dict=None, length_dep=True, gap_penalty=-10):
    X = construct_distance_vector(Xlabels, length_dep, gap_penalty)
    Z = linkage(X, 'average')

    cluster_arr = fcluster(Z, max_distance, criterion='distance')

    clusters = clusters_to_dict(cluster_arr, Xlabels)
    return clusters


def cluster_data_and_generate_dendrogram(Xlabels, max_distance=0.2, colors_dict=False, structural_cluster_dict=None, length_dep=True, gap_penalty=-10):
    X = construct_distance_vector(Xlabels, length_dep, gap_penalty)
    Z = linkage(X, 'average')

    for i in range(len(Xlabels)):
        print(str(i) + ": " + Xlabels[i])

    plt.figure(figsize=(10, 5))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')

    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        labels=Xlabels
    )

    if colors_dict:
        num_to_color = {1: 'r', 2: 'g', 3: 'b', 4: 'm', 5: 'y', 6: 'c'}
        ax = plt.gca()

        x_labels = ax.get_xticklabels()
        for x in x_labels:
            if structural_cluster_dict[x.get_text()] != 0:
                x.set_color(num_to_color[structural_cluster_dict[x.get_text()]])

    plt.show()

    cluster_arr = fcluster(Z, max_distance, criterion='distance')

    clusters = clusters_to_dict(cluster_arr, Xlabels)
    return clusters

def create_motif_from_fasta_file(fasta_filename, out_filename, generate_pssm=False):
    instances = []
    with open(fasta_filename) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            instances.append(Seq(seq, IUPAC.protein))

    m = motifs.create(instances, IUPAC.protein)
    m.weblogo(out_filename, format='svg')

    if generate_pssm:
        pssm_file = open(out_filename[:-10] + "_pssm.txt", "w+")
        for i in range(len(m.pwm)):
            for j in range(len(m.pwm[i])):
                pssm_file.write(str(m.pwm[i][j]))
                pssm_file.write("\t")

            pssm_file.write("\n")
        pssm_file.close()


# takes in the clusters
# generates the correctly formatted text file of the CDRs in the cluster
# if a size is given, only clusters of that size will be written to files
def generate_motifs(clusters, file_path, clust_size=None):
    for key in clusters:
        print ("processing cluster " + str(key))
        process_cluster = True
        if clust_size:
            size = len(clusters[key])
            if size != clust_size:
                process_cluster = False
        if process_cluster:
            file_name = file_path + str(key) + ".txt"
            f = open(file_name, "w")
            clust = clusters[key]
            i = 0
            for read in clust:
                f.write(">" + str(i) + "\n")
                f.write(read + "\n\n")
                i += 1
            f.close()

            create_motif_from_fasta_file(file_name, file_name[:-4] + "_motif.svg")

def generate_motif_from_cluster(key, cluster, file_path, clust_size=None):
    clusters = {key: cluster}
    generate_motifs(clusters, file_path, clust_size)

# n = size of clusters in smallest distance clustering, to keep track of
# clusters2 = cluster_data(Xlabels, max_distance=0.2)
# clusters3 = cluster_data(Xlabels, max_distance=0.3)
# clusters4 = cluster_data(Xlabels, max_distance=0.4)
# default is to find the cluster progression of all of the clusters in clusters2
# can also pass in a cluster in clust, which will return the progression of a specific cluster
def cluster_progression(out_dir, clusters2, clusters3, clusters4, n=None, clust=None):

    motif_names = []

    if clust != None:
        for c in clusters2:
            if clust == clusters2[c]:
                 # process cluster - generate motif
                motif_name = str(c) + "_0.2"
                generate_motif_from_cluster(motif_name, clusters2[c], out_dir)
                motif_names.append(out_dir + motif_name + "_pssm.txt")


                # get first CDR - all unique, so can use to find how the cluster evolved
                first_cdr = clusters2[c][0]

                for j in clusters3:
                    if first_cdr in clusters3[j]:
                        motif_name = str(c) + "_0.3"
                        generate_motif_from_cluster(motif_name, clusters3[j], out_dir)
                        motif_names.append(out_dir + motif_name + "_pssm.txt")

                        break

                for i in clusters4:
                    if first_cdr in clusters4[i]:
                        motif_name = str(c) + "_0.4"
                        generate_motif_from_cluster(motif_name, clusters4[i], out_dir)
                        motif_names.append(out_dir + motif_name + "_pssm.txt")

                        break

    else:
        f = open(out_dir + "cluster_progression_results.txt", "w")
        correct_n = True
        # extract all clusters of correct size and process
        for c in clusters2:
            if n:
                if len(clusters2[c]) == n:
                    correct_n = True
                else:
                    correct_n = False
            if correct_n and len(clusters2[c]) >= 2:
                f.write("\n\noriginal cluster, at distance 0.2: ")
                f.write(str(clusters2[c]))
                f.write("\ncluster number " + str(c))

                # process cluster - generate motif
                motif_name = str(c) + "_0.2"
                generate_motif_from_cluster(motif_name, clusters2[c], out_dir)
                motif_names.append(out_dir + motif_name + "_pssm.txt")

                # get first CDR - all unique, so can use to find how the cluster evolved
                first_cdr = clusters2[c][0]

                for j in clusters3:
                    if first_cdr in clusters3[j]:
                        f.write("\ncluster at distance 0.3: ")
                        f.write(str(clusters3[j]))

                        motif_name = str(c) + "_0.3"
                        generate_motif_from_cluster(motif_name, clusters3[j], out_dir)
                        motif_names.append(out_dir + motif_name + "_pssm.txt")

                        break

                for i in clusters4:
                    if first_cdr in clusters4[i]:
                        f.write("\ncluster at distance 0.4: ")
                        f.write(str(clusters4[i]))

                        motif_name = str(c) + "_0.4"
                        generate_motif_from_cluster(motif_name, clusters4[i], out_dir)
                        motif_names.append(out_dir + motif_name + "_pssm.txt")

                        break
        f.close()
    return motif_names

