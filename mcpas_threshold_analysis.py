# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from nearest_neighbor_threshold_analysis import run_threshold_analysis_no_plot
import csv, sys

if len(sys.argv) < 4:
    print("need the following arguments: mcpas data file path, output_directory, "
          "small results filename, large results filename")
    sys.exit()

# data file acquired at http://friedmanlab.weizmann.ac.il/McPAS-TCR/ -> download the complete database
mcpas_data_file = sys.argv[1] # datasets/mcpas_data/McPAS-TCR.csv
output_dir = sys.argv[2]
if output_dir[-1] != '/': output_dir += '/'
large_file_name = str(sys.argv[3]) # large
if large_file_name[-1] != '_': large_file_name += '_'
small_file_name = sys.argv[4] # small
if small_file_name[-1] != '_': small_file_name += "_"

mcpas_data = open(mcpas_data_file, "rb")
mcpas_reader = csv.reader(mcpas_data)

CDR3b_i = 1
Species_i = 2
Pathology_i = 4

human_pathology_dict = {}

abbrevs = {"Multiple sclerosis (MS)" : "MS",
           "Clear cell renal carcinoma": "clear_cell",
           "Yellow fever virus" : "yellow_fever",
            "Colorectal cancer" : "colorectal",
            "Rheumatoid Arthritis (RA)": "RA",
            "Celiac disease" : "celiac",
            "Hepatitis C virus" : "HepC",
            "Herpes simplex virus 2 (HSV2)" : "HSV2",
            "Diabetes Type 1" : "diabetes",
            "Human immunodeficiency virus (HIV)" : "HIV",
            "Epstein Barr virus (EBV)" : "EBV",
            "Cytomegalovirus (CMV)" : "CMV"
           }

for row in mcpas_reader:
    # only want to analyze human
    if row[Species_i] == "Human":
        pathology = row[Pathology_i]
        CDR3b = row[CDR3b_i]

        if pathology in human_pathology_dict:
            human_pathology_dict[pathology].append(CDR3b[1:])
        else:
            human_pathology_dict[pathology] = [CDR3b[1:]]

small_paths = {}
large_paths = {}

for pathology in human_pathology_dict:

    seqs = list(set(human_pathology_dict[pathology]))   # remove duplicates
    size = len(seqs)

    print(pathology)
    print("size: " + str(size))

    if size < 50:
        print("too small, don't consider")
    else:
        if pathology in abbrevs:
            name = abbrevs[pathology]
        else:
            name = pathology
        if size > 50 and size < 300:
            print("small repertoire")
            small_paths[name] = seqs
        else:
            print("large repertoire")
            large_paths[name] = seqs
    print('')


def set_up_results_csv(out_path):
    filename = out_path + "results.csv"
    with open(filename, mode='w') as outfile:
        csv_writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['p1', 'p2', 'threshold', 'unidentified', 'correct', 'incorrect', 'unidentified_frac', 'correct_frac', 'incorrect_frac'])

    return filename

def get_mcpas_results(repertoire1_name, repertoire2_name, pathology_dict, output_file):
    results = run_threshold_analysis_no_plot([repertoire1_name, repertoire2_name], pathology_dict)

    r0 = results[0]
    r2 = results[0.2]
    r3 = results[0.3]
    r4 = results[0.4]
    r1 = results[1]

    count = float(r0[0])

    with open(output_file, mode='a') as outfile:
        csv_writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        csv_writer.writerow([repertoire1_name, repertoire2_name, 0, r0[0], r0[1], r0[2], r0[0]/count, r0[1]/count, r0[2]/count])
        csv_writer.writerow([repertoire1_name, repertoire2_name, 0.2, r2[0], r2[1], r2[2], r2[0]/count, r2[1]/count, r2[2]/count])
        csv_writer.writerow([repertoire1_name, repertoire2_name, 0.3, r3[0], r3[1], r3[2], r3[0]/count, r3[1]/count, r3[2]/count])
        csv_writer.writerow([repertoire1_name, repertoire2_name, 0.4, r4[0], r4[1], r4[2], r4[0]/count, r4[1]/count, r4[2]/count])
        csv_writer.writerow([repertoire1_name, repertoire2_name, 1, r1[0], r1[1], r1[2], r1[0]/count, r1[1]/count, r1[2]/count])


output_file_small = set_up_results_csv(output_dir + small_file_name)
small_paths_list = small_paths.keys()

for i in range(len(small_paths_list)):
    for j in range(i + 1, len(small_paths_list)):
        get_mcpas_results(small_paths_list[i], small_paths_list[j], small_paths, output_file_small)


output_file_large = set_up_results_csv(output_dir + large_file_name)
large_paths_list = large_paths.keys()

for i in range(len(large_paths_list)):
    for j in range(i + 1, len(large_paths_list)):
        get_mcpas_results(large_paths_list[i], large_paths_list[j], large_paths, output_file_large)

