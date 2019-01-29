from __future__ import print_function
import csv, collections, os, sys

# args
# 1: cdr -- a3 or b3
# 2: stcr_data directory -- ...stcr_data/ -- should contain <cdr>/details/*.tsv
# 3: IMGT-renumbered structures directory -- ...all_structures/imgt/
# ex: python get_cdrs_stcr.py b3 ../datasets/stcr_data/ ~/Downloads/all_structures/imgt/

cdr = sys.argv[1]
root = sys.argv[2]
if root[-1] != '/': root += '/'
structs = sys.argv[3]
if structs[-1] != '/': structs += '/'

if not os.path.exists(root+cdr+'/structs'): os.mkdir(root+cdr+'/structs')

# http://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
# cdr:(start,end) as in range()
numbering = {'b3':(105,118), 'a3':(105,118)}

# convert 3-char to 1-char code
AA31 = dict([('ALA','A'),('CYS','C'),('ASP','D'),('GLU','E'),('PHE','F'),('GLY','G'),('HIS','H'),('ILE','I'),('LYS','K'),('LEU','L'),('MET','M'),('ASN','N'),('PRO','P'),('GLN','Q'),('ARG','R'),('SER','S'),('THR','T'),('VAL','V'),('TRP','W'),('TYR','Y')])

# all the canonical clusters for this cdr, in tsv files downloaded from stcrdab
groups = [g[:-4] for g in os.listdir(root+cdr+'/details') if g[-4:]=='.tsv']

# which residues to extract for this cdr
resi_range = set(range(*numbering[cdr]))

cdrs = [] # list of cdr specs to be saved out
pdb_seqs = collections.defaultdict(set) # pdb_id => cdr sequences, for duplicate checking
for group in groups:
    for record in csv.DictReader(open(root+cdr+'/details/'+group+'.tsv'), delimiter='\t'):
        pdb_id = record['pdb']
        pdb_chain = record[cdr[0].upper()+'chain']
        if pdb_chain == 'NA': continue # no structure?
        out_pdb_id = pdb_id+pdb_chain+'_'+cdr
        atoms = [] # list of ATOM records within the cdr
        aas = {} # AA sequence of the cdr, by position
        # extract atoms and seq
        for row in open(structs+pdb_id+'.pdb'):
            if row[:4]=='ATOM' and row[21]==pdb_chain and int(row[23:26]) in resi_range:
                atoms.append(row)
                aas[int(row[23:26])] = AA31[row[17:20]]
        seq = ''.join(aas[pos] for pos in sorted(aas))
        pdb_seqs[pdb_id].add(seq)
        cdrs.append((group,len(seq),seq,pdb_id,pdb_chain))
        # save pdb file of just the cdr
        with open(root+cdr+'/structs/'+out_pdb_id+'.pdb','w') as loopfile:
            for atom in atoms:
                loopfile.write(atom)

# note when different cdrs (presumably some missing electron density)
with open(root+cdr+'/log.txt','w') as logfile:
    for pdb_id,seqs in pdb_seqs.items():
        if len(seqs)>1:
            print('different cdrs in',pdb_id,seqs, file=logfile)

# save all the cdr specs
with open(root+cdr+'/cdrs.csv','w') as outfile:
    outcsv = csv.writer(outfile)
    outcsv.writerow(('group','len','seq','pdb','chain'))
    outcsv.writerows(cdrs)
