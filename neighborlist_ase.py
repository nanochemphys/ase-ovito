from ase.io import read, write
from ase.calculators.emt import EMT
from ase import Atom, Atoms
from ase import neighborlist
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.data import covalent_radii

## create ase-neighbor list, natural cuttoffs are multiple covalent radii
atoms = read('clusters.xyz')
natoms = atoms.get_global_number_of_atoms()
cutoff = neighborlist.natural_cutoffs(atoms, 2.3)
nl = NeighborList(cutoff, skin=0, self_interaction=True, bothways=True)
nl.update(atoms)
coord = []
for i, atom in enumerate(atoms):
	temp_clust = []
	x0 = atoms.positions[i][0]
	y0 = atoms.positions[i][1]
	z0 = atoms.positions[i][2]
	indices, offsets = nl.get_neighbors(i)
	for j, offset in zip(indices, offsets):
		temp_clust.append(atoms.positions[j])
	temp_clust = [temp for temp in temp_clust if (temp[0] != x0 and temp[1] != y0 and temp[2] != z0)]
	coord.append(len(temp_clust))
## write to out file
file = open('cluster_nl.xyz', 'w')
file1.write('%s\n\n' % natoms)
for i in range(0, natoms, 1):
	if int(tags[i]) == 1:
		file1.write('%s\t' % atoms.symbols[i])
		file1.write('%s\t'  % atoms.positions[i][0])
		file1.write('%s\t'  % atoms.positions[i][1])
		file1.write('%s\t'  % atoms.positions[i][2])
		file1.write('%s\n' % coord[i])
file1.close()


    
