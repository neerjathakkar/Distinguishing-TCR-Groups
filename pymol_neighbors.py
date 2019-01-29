# By Chris Bailey-Kellogg for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

import csv, subprocess, sys

cdr_id = sys.argv[1] # b3 or a3
target = sys.argv[2] # cdr sequence
snbr = sys.argv[3]   # sequence of closest neighbor in its cluster
onbr = sys.argv[4]   # sequence of closest neighbor not in its cluster

base = 'datasets/stcr_data/' # could allow this to be an arg...

# load mapping from sequences to representative structures
cdr2struct = {}
for row in csv.DictReader(open(base+cdr_id+'_reps.csv')):
    cdr2struct[row['cdr']] = row['file']

if target not in cdr2struct: raise Exception('no struct for target '+target)
if snbr not in cdr2struct: raise Exception('no struct for snbr '+snbr)
if onbr not in cdr2struct: raise Exception('no struct for onbr '+onbr)

# pymol commands
script =  """load {structs}/{target}, target
load {structs}/{snbr}, snbr
load {structs}/{onbr}, onbr
align snbr, target
align onbr, target
hide all; show cartoon
cartoon tube
set cartoon tube_radius 0.3
color blue, target
color green, snbr
color orange, onbr
bg_color white
set opaque_background, 0
set depth_cue, 0
reset""" .format(structs=base+cdr_id+'/structs', target=cdr2struct[target], snbr=cdr2struct[snbr], onbr=cdr2struct[onbr])

# run pymol
sp = subprocess.Popen('pymol -p -K', stdin=subprocess.PIPE, shell=True)
(sp_stdout, sp_stderr) = sp.communicate(script)
