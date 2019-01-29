# By Chris Bailey-Kellogg for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

import csv, collections, sys
from superimpose import load_backbone, superimposed_rmsd

# args:
# 1: cdr id -- a3 or b3
# [2]: stcr_data directory (defaults to datasets/stcr_data/)
# [3]: stcr_results directory (defaults to results/stcr_results/)

cdr_id = sys.argv[1]
if len(sys.argv)>2:
    data_dir = sys.argv[2]
    if data_dir[-1] != '/': data_dir += '/'
else:
    data_dir = 'datasets/stcr_data/'
if len(sys.argv)>3:
    results_dir = sys.argv[3]
    if results_dir[-1] != '/': results_dir += '/'
else:
    results_dir = 'results/stcr_results/'

cdr_files = collections.defaultdict(list)
for row in csv.DictReader(open(data_dir+cdr_id+'/cdrs.csv')):
    cdr_files[row['seq']].append(row['pdb']+row['chain']+'_'+cdr_id+'.pdb')

# when multiple struct files, choose one with smallest total rmsd to others
rep_structs = {}
rep_details = []
for cdr,files in cdr_files.iteritems():
    if len(files)==1:
        rep_structs[cdr] = load_backbone(data_dir+cdr_id+'/structs/'+files[0])
        rep_details.append((cdr,files[0], 1, 0))
    else:
        structs = [load_backbone(data_dir+cdr_id+'/structs/'+f) for f in files]
        dists = [0]*len(structs)
        best = None; best_dist = None
        for i in range(len(structs)):
            for j in range(i+1, len(structs)):
                d = superimposed_rmsd(structs[i], structs[j])
                dists[i] += d
                dists[j] += d
            if best is None or dists[i] < best_dist:
                best = i; best_dist = dists[i]
        rep_structs[cdr] = structs[i]
        rep_details.append((cdr,files[i], len(structs), best_dist/len(structs)))
with open(data_dir+cdr_id+'_reps.csv', 'w') as outfile:
    outcsv = csv.writer(outfile)
    outcsv.writerow(['cdr','file','n','avg rmsd'])
    outcsv.writerows(rep_details)

# rmsds for specified pairs
rmsds = []
missings = set()
for row in csv.DictReader(open(results_dir+cdr_id+'_stcr_pairs.csv')):
    if row['cdr'] not in rep_structs:
        missings.add(row['cdr'])
        continue
    if row['nbr'] not in rep_structs:
        missings.add(row['nbr'])
        continue
    rmsds.append((row['cdr'],row['nbr'], superimposed_rmsd(rep_structs[row['cdr']], rep_structs[row['nbr']])))
with open(results_dir+cdr_id+'_pair_rmsds.csv', 'w') as outfile:
    outcsv = csv.writer(outfile)
    outcsv.writerow(['cdr','nbr','rmsd'])
    outcsv.writerows(rmsds)
with open(results_dir+cdr_id+'_missings.txt', 'w') as outfile:
    for cdr in missings:
        outfile.write(cdr+'\n')
