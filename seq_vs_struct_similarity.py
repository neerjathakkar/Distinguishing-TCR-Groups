# By Chris Bailey-Kellogg for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

import csv, sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 

# args:
# 1: cdr id -- a3 or b3
# [2]: stcr_results directory (defaults to results/stcr_results/)

cdr_id = sys.argv[1]
if len(sys.argv)>2:
    results_dir = sys.argv[2]
    if results_dir[-1] != '/': results_dir += '/'
else:
    results_dir = 'results/stcr_results/'

# (cdr, nbr) -> distance, for both sequence and structural distances
seq_dists = {}
struct_dists = {}
# seq -> nearest seq nbr in same or other cluster
same = {}
other = {}
for row in csv.DictReader(open(results_dir+cdr_id+'_stcr_pairs.csv')):
    seq_dists[row['cdr'],row['nbr']] = float(row['seq dist'])
    if row['cluster'] == row['nbr cluster']: same[row['cdr']] = row['nbr']
    else: other[row['cdr']] = row['nbr']
for row in csv.DictReader(open(results_dir+cdr_id+'_pair_rmsds.csv')):
    struct_dists[row['cdr'],row['nbr']] = float(row['rmsd'])

assert(set(same.keys()) == set(other.keys()))

# save out compiled distances
with open(results_dir+cdr_id+'_seq_vs_struct.csv','w') as outfile:
    csvout = csv.writer(outfile)
    csvout.writerow(['cdr','correct','snbr','snbr seq','snbr struct','onbr','onbr seq','onbr struct','delta seq','delta struct'])
    for cdr in same:
        snbr = same[cdr]; onbr = other[cdr]
        csvout.writerow([cdr, 1 if seq_dists[cdr,snbr]<=seq_dists[cdr,onbr] else 0,
                             snbr, seq_dists[cdr,snbr], struct_dists.get((cdr,snbr),''),
                             onbr, seq_dists[cdr,onbr], struct_dists.get((cdr,onbr),''),
                             seq_dists[cdr,snbr] - seq_dists[cdr,onbr],
                             struct_dists[cdr,snbr] - struct_dists[cdr,onbr] if ((cdr,snbr) in struct_dists and (cdr,onbr) in struct_dists) else ''])

# for each cdr, plot delta seq (same vs. other nbr) vs. delta struct
# different color for correct vs. incorrect classifications, though that's redundant since quadrant makes that clear
cmap = plt.cm.get_cmap('gnuplot2') # see https://matplotlib.org/examples/color/colormaps_reference.html
fig = plt.figure()
plt.axvline(x=0, color='gray', linestyle='--')
plt.axhline(y=0, color='gray', linestyle='--')
for cdr in same:
    snbr = same[cdr]; onbr = other[cdr]
    if (cdr,snbr) not in struct_dists:
        print 'no struct_dist for',cdr,snbr
        continue
    if (cdr,onbr) not in struct_dists:
        print 'no struct_dist for',cdr,onbr
        continue
    seq_delta = seq_dists[cdr,snbr]-seq_dists[cdr,onbr]
    struct_delta = struct_dists[cdr,snbr]-struct_dists[cdr,onbr]
    plt.plot(seq_delta, struct_delta, '.', markersize=12, color=cmap(seq_dists[cdr,snbr]))
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('$\Delta$ seq', fontsize=18)
plt.ylabel('$\Delta$ struct', fontsize=18)
divider = make_axes_locatable(plt.gca())
cax = divider.new_horizontal(size="5%", pad=0.2, pack_start=False) # thanks to https://stackoverflow.com/questions/13310594/positioning-the-colorbar
fig.add_axes(cax)
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical')
cb.set_label('seq dist', fontsize=15)
plt.tight_layout()
plt.savefig(results_dir+cdr_id+'_seq_vs_struct.pdf')
