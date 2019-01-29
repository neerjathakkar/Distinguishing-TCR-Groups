# By Chris Bailey-Kellogg for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

import numpy as np
import sys

def load_backbone(filename):
  """Loads just the backbone atom records from the PDB file. Minimal error checking!"""
  atoms = []
  for line in open(filename, 'r'):
	if line[0:4] == 'ATOM' and (line[13:15]=='N ' or line[13:15]=='CA' or line[13:15]=='C '):
	  atoms.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
  return atoms

def superimpose_onto(mobile_atoms, fixed_atoms):
  """Finds a transformation for the mobile atoms onto the fixed atoms in order to minimize RMSD between corresponding pairs (lists must be of the same length).
  Returns a function to do the transformation."""
  t1 = np.average(mobile_atoms, 0)
  t2 = np.average(fixed_atoms, 0)
  centered_pairs = [(a1-t1, a2-t2) for (a1,a2) in zip(mobile_atoms, fixed_atoms)]

  C = np.zeros((3,3))
  for d1 in xrange(3):
	for d2 in xrange(3):
	  C[d2][d1] = np.sum(a1[d1]*a2[d2] for (a1,a2) in centered_pairs)

  U,W,VT = np.linalg.svd(C)
  D = np.identity(3); D[2][2] = np.linalg.det(np.dot(U,VT))
  R = np.dot(U, np.dot(D,VT))

  return lambda a: np.dot(R, a-t1) + t2

def transform_pdb(xform, in_filename, out_filename):
  """Produces in out_filename a copy of the ATOMs and HETATMs in in_filename, transformed according to xform. Echoes other lines."""
  with open(out_filename, 'w') as out:
	for line in open(in_filename,'r'):
	  if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
		p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
		pp = xform(p)
		out.write(line[:30] + str(round(pp[0],3)).rjust(8) + str(round(pp[1],3)).rjust(8) + str(round(pp[2],3)).rjust(8) + line[54:])
	  else:
		out.write(line)

def rmsd(atoms1, atoms2):
  """RMSD between corresponding atoms in the two lists (must be of the same length)"""
  return np.sqrt(np.average([np.square(np.linalg.norm(a1-a2)) for (a1,a2) in zip(atoms1,atoms2)]))

def superimposed_rmsd(atoms1, atoms2):
  """RMSD between corresponding atoms after finding the berst superimposition"""
  xform = superimpose_onto(atoms1, atoms2)
  return rmsd([xform(a) for a in atoms1], atoms2)

if __name__=='__main__':
  if len(sys.argv) != 4:
	print 'usage: <mobile pdb in> <fixed pdb> <mobile pdb out>'
	exit()
  mobile_atoms = load_backbone(sys.argv[1])
  fixed_atoms = load_backbone(sys.argv[2])
  if len(mobile_atoms) != len(fixed_atoms):
	print 'different lengths', len(mobile_atoms), len(fixed_atoms)
  xform = superimpose_onto(mobile_atoms, fixed_atoms)
  print 'ca rmsd', rmsd([xform(a) for a in mobile_atoms], fixed_atoms)
  transform_pdb(xform, sys.argv[1], sys.argv[3])
