# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from nearest_neighbor_threshold_analysis import run_threshold_analysis_no_plot
import csv, sys

if len(sys.argv) < 3:
    print("need the following arguments: glanville data file path, output_directory, "
          "results filename")
    sys.exit()

mcpas_data_file = sys.argv[1] # datasets/glanville_data/glanville_s2_curated_data.csv
output_dir = sys.argv[2]
if output_dir[-1] != '/': output_dir += '/'
output_file_name = str(sys.argv[3])
if output_file_name[-1] != '_': output_file_name += '_'

glanville_data = open(mcpas_data_file, "rb")
glanville_reader = csv.reader(glanville_data)

CDR3b_i = 4
antigen_i = 5
hla_i = 2

antigen_dict = {}

i = 0

for row in glanville_reader:
    antigen = row[antigen_i]
    if antigen != 'Antigen' and antigen != '':
        if antigen == "pp65":
            if row[hla_i] == "HLA-A2":
                antigen = "pp65_A2"
            elif row[hla_i] == "HLA-B7":
                antigen = "pp65_B7"

        CDR3b = row[CDR3b_i]

        if antigen in antigen_dict:
            antigen_dict[antigen].append(CDR3b[1:])
        else:
            antigen_dict[antigen] = [CDR3b[1:]]

antigens = {}

for antigen in antigen_dict:
    print(antigen)
    antigens[antigen] = list(set(antigen_dict[antigen])) # remove duplicates

    size = len(antigens[antigen])
    print("size: " + str(size))

# print(antigens)

def set_up_results_csv(out_path):
    filename = out_path + "results.csv"
    with open(filename, mode='w') as outfile:
        csv_writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['p1', 'p2', 'threshold', 'unidentified', 'correct', 'incorrect', 'unidentified_frac', 'correct_frac', 'incorrect_frac'])

    return filename

def get_antigen_results(repertoire1_name, repertoire2_name, pathology_dict, output_file):
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


output_file_small = set_up_results_csv(output_dir + output_file_name)
antigens_list = antigens.keys()

for i in range(len(antigens_list)):
    for j in range(i + 1, len(antigens_list)):
        get_antigen_results(antigens_list[i], antigens_list[j], antigens, output_file_small)


