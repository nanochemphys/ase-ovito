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
##Translation from DL POLY-format to xyz-format
with open('REVCON', 'r') as file:
    strngs = file.readlines()
file.close()
nstr = len(strngs)
tmp = strngs[1].split()
natoms = int(tmp[2])
tmp = strngs[2].split()
cell_x = float(tmp[0])
tmp = strngs[3].split()
cell_y = float(tmp[1])
tmp = strngs[4].split()
cell_z = float(tmp[2])
atom_symb = []
atom_x = []
atom_y = []
atom_z = []
velc_x = []
velc_y = []
velc_z = []
forc_x = []
forc_y = []
forc_z = []
for i in range (5, nstr, 4):
    tmp = strngs[i].split()
    atom_symb.append(tmp[0])
    tmp = strngs[i+1].split()
    atom_x.append(tmp[0])
    atom_y.append(tmp[1])
    atom_z.append(tmp[2])
    tmp = strngs[i+2].split()
    velc_x.append(tmp[0])
    velc_y.append(tmp[1])
    velc_z.append(tmp[2])
    tmp = strngs[i+3].split()
    forc_x.append(tmp[0])
    forc_y.append(tmp[1])
    forc_z.append(tmp[2])
with open('cluster_tmp.xyz', 'w') as file:
    file.write(str(natoms) + '\n\n')
    for i in range(natoms):
        file.write(atom_symb[i] + '\t' + atom_x[i] + '\t' + atom_y[i] + '\t' + atom_z[i] + '\n')
file.close()

## Create surface layer by OVITO-tools
pipeline = import_file('cluster_tmp.xyz')
pipeline.modifiers.append(ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 3.1,
    select_surface_particles = True))
export_file(pipeline, 'clusters_out.xyz', 'xyz', columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "Selection"])
pipeline.modifiers.append(ClusterAnalysisModifier(
    cutoff=3.1, 
    compute_gyration = True))
data = pipeline.compute()
cluster_table = data.tables['clusters']
with open("data_cluster.txt", 'w') as file:
    file.write(str(cluster_table['Radius of Gyration'][...]))

##Read xyz-file with ovito-tags 
atoms = read('clusters_out.xyz')
at_symb = atoms.symbols[0]
cell0 = atoms.cell
calc = EMT()
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
## create ase-neighbor list, natural cuttoffs are multiple covalent radii
cutoff = neighborlist.natural_cutoffs(atoms, 2.3)
nl = NeighborList(cutoff, skin=0, self_interaction=True, bothways=True)
nl.update(atoms)
coord = []
ro_fs = []
ro_sc = []
ro_gu = []
## Calculate electron density by Gupta, Sutton-Chen, Finnis-Sinclair potential
for i, atom in enumerate(atoms):
	temp_clust = []
	indices, offsets = nl.get_neighbors(i)
	for j, offset in zip(indices, offsets):
		temp_clust.append(atoms.positions[j])
	coord.append(len(temp_clust))
	x0 = atoms.positions[i][0]
	y0 = atoms.positions[i][1]
	z0 = atoms.positions[i][2]
	temp_clust = [temp for temp in temp_clust if (temp[0] != x0 and temp[1] != y0 and temp[2] != z0)]
	ro_fs_sum = 0
	ro_sc_sum = 0
	ro_gu_sum = 0
	for l, temp in enumerate(temp_clust):
		rij = sqrt((x0-temp[0])**2 + (y0-temp[1])**2 + (z0-temp[2])**2)
		ro_sc_i = (a/rij)**m
		ro_sc_sum = ro_sc_sum + ro_sc_i
		ro_gu_i = exp(-2*q*((rij-r0)/r0))
		ro_gu_sum = ro_gu_sum + ro_gu_i
        if r_ij-d <= d:
            ro_fs_i = (A**2)*((rij-d)**2 + (B**2)*(rij-d)**4)
            ro_fs_sum = ro_fs_sum + ro_fs_i
        else:
            ro_fs_i = 0
            ro_fs_sum = ro_fs_sum + ro_fs_i
	ro_fs.append(ro_fs_sum)
	ro_sc.append(ro_sc_sum)
	ro_gu.append(ro_gu_sum)
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
		file1.write('%s\t' % ro_fs[i])
		file1.write('%s\t' % ro_sc[i])
		file1.write('%s\n' % ro_gu[i])
file1.close()


    
