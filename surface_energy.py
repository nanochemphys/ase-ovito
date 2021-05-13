from ase.io import read, write
from ase.calculators.emt import EMT
from ase import Atom, Atoms
from ase import neighborlist
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.data import covalent_radii
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier
from ovito.modifiers import ConstructSurfaceModifier
from math import sqrt, exp
import numpy as np

np.seterr(divide='ignore', invalid='ignore')
## Semiempirical coefficients Au
d = 4.46
A = 0.0137025
B = -53.9630
a = 4.0800
m = 8
r0 = 2.8838
q = 3.9472
## Create surface layer by OVITO-tools
pipeline = import_file('cluster_tmp.xyz')
pipeline.modifiers.append(ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 3.1,
    select_surface_particles = True))
export_file(pipeline, 'clusters_out.xyz', 'xyz', columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "Selection"])
##Read xyz-file with ovito-tags 
atoms = read('clusters_out.xyz')
cell0 = atoms.cell
calc = EMT()
at_symb = atoms.symbols[0]
f = open('clusters_out.xyz', 'r')
lines = f.readlines()
while len(lines) > 0:
	tags = []
	surf_nn = 0
	natoms = int(lines.pop(0))
	lines.pop(0)
	for _ in range(natoms):
		line = lines.pop(0)
		symbol, x, y, z, tag = line.split()[:5]
		surf_nn = surf_nn + int(tag)
		tags.append(int(tag))
## create ase-neighborlist, natural cuttoffs are multiple covalent radii
cutoff = neighborlist.natural_cutoffs(atoms, 2.3)
nl = NeighborList(cutoff, skin=0, self_interaction=True, bothways=True)
nl.update(atoms)
surf_energy = []
surf_energy_weight = []
coord = []
e_il_sum = 0
## Calculate surface energy
for i, atom in enumerate(atoms):
	temp_clust = []
	indices, offsets = nl.get_neighbors(i)
	for j, offset in zip(indices, offsets):
		temp_clust.append(atoms.positions[j]) ## append to list
	atoms_i = atoms.positions[i]
	temp_clust = [temp for temp in temp_clust if (temp[0] != atoms_i[0] and temp[1] != atoms_i[1] and temp[2] != atoms_i[2])]
	for l, temp in enumerate(temp_clust):
		pair_atoms = Atoms('Au2', positions=[(atoms_i[0], atoms_i[1], atoms_i[2]),(temp[0],temp[1], temp[2])], cell=cell0, pbc=[0,0,0])
		pair_atoms.set_calculator(calc) 
		e_il = pair_atoms.get_potential_energy()
		e_il_sum = e_il_sum + e_il
	e_surf = e_il_sum/(len(temp_clust))
	surf_energy.append(e_surf)
	coord.append(len(temp_clust))
k_max = max(coord)
for i, atom in enumerate(atoms):
	e_surf_weight = k_max * surf_energy[i] / coord[i]
	surf_energy_weight.append(e_surf_weight)
## write to out file
file1 = open('%s_clust%s_surf%s.xyz' % (at_symb, natoms, surf_nn), 'w')
file1.write('%s\n\n' % surf_nn)
for i in range(0, natoms, 1):
	if int(tags[i]) == 1:
		file1.write('%s\t' % atoms.symbols[i])
		file1.write('%s\t'  % atoms.positions[i][0])
		file1.write('%s\t'  % atoms.positions[i][1])
		file1.write('%s\t'  % atoms.positions[i][2])
		file1.write('%s\t' % coord[i])
		file1.write('%s\t' % surf_energy[i])
		file1.write('%s\n' % surf_energy_weight[i])
file1.close()


    
