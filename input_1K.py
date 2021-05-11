from ase.io import read, write
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier
from ovito.modifiers import ConstructSurfaceModifier
import numpy as np
from ase import Atom, Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, Hookean
from ase.visualize import view
from ase import neighborlist
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.data import covalent_radii
from math import sqrt, exp

with open('REVCON_1K', 'r') as file:
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
with open('cluster_tmp_1K.xyz', 'w') as file:
    file.write(str(natoms) + '\n\n')
    for i in range(natoms):
        file.write(atom_symb[i] + '\t' + atom_x[i] + '\t' + atom_y[i] + '\t' + atom_z[i] + '\n')
file.close()

pipeline = import_file('cluster_tmp_1K.xyz')

pipeline.modifiers.append(ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 3.1,
    select_surface_particles = True))

# Export results of the clustering algorithm to a text file:
export_file(pipeline, 'clusters_out_1K.xyz', 'xyz', columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "Selection"])

pipeline.modifiers.append(ClusterAnalysisModifier(
    cutoff=3.1, 
    compute_gyration = True))

# Directly access information stored in the DataTable:
data = pipeline.compute()
cluster_table = data.tables['clusters']
with open("data_cluster_1K.txt", 'w') as file:
    file.write(str(cluster_table['Radius of Gyration'][...]))

##Read cluster
atoms = read('clusters_out_1K.xyz')
f = open('clusters_out_1K.xyz', 'r')
lines = f.readlines()
while len(lines) > 0:
   tags = []
   surf_nn = 0
   natoms = int(lines.pop(0))
   lines.pop(0)  # Comment line; ignored
   for _ in range(natoms):
       line = lines.pop(0)
       symbol, x, y, z, tag = line.split()[:5]
       surf_nn = surf_nn + int(tag)
       tags.append(int(tag))
calc = EMT()
nn = atoms.get_global_number_of_atoms()
cell0 = atoms.cell
cutoff = neighborlist.natural_cutoffs(atoms, 1.2)
## create neighbor list, skin = 0.0 - 0.2
nl = NeighborList(cutoff, skin=0, self_interaction=True, bothways=True)
nl.update(atoms)
d = 4.46
A = 0.0137025
B = -53.9630
a = 4.0800
m = 8
r0 = 2.8838
q = 3.9472
## arrays for list of neighbor atoms positions and extra surface energy
surf_energy = []
coord = []
ro_fs = []
ro_sc = []
ro_gu = []
##create list 
for i in range(0, nn, 1):
   temp_clust = []
   indices, offsets = nl.get_neighbors(i)
   for j, offset in zip(indices, offsets):
       temp_clust.append(atoms.positions[j]) ## append to list
   k = len(temp_clust)-1
#   print(k)
   e_il_sum = 0
## Calculate energy of pair atoms
   for l in range(0, k, 1):
       pair_atoms = Atoms('Au2', positions=[(atoms.positions[i][0], atoms.positions[i][1], atoms.positions[i][2]),(temp_clust[l][0],temp_clust[l][1], temp_clust[l][2])], cell=cell0, pbc=[0,0,0])
       pair_atoms.set_calculator(calc) 
       e_il = pair_atoms.get_potential_energy()
       e_il_sum = e_il_sum + e_il
   e_surf = e_il_sum/k
   surf_energy.append(e_surf)
   coord.append(k)
k_max = max(coord)
surf_energy_weight = []
for i in range(0, nn, 1):
   e_surf_weight = k_max * surf_energy[i] / coord[i]
   surf_energy_weight.append(e_surf_weight)
for i in range(0, nn, 1):
   temp_clust = []
   indices, offsets = nl.get_neighbors(i)
   for j, offset in zip(indices, offsets):
       temp_clust.append(atoms.positions[j]) ## append to list
   h = len(temp_clust)
   x0 = atoms.positions[i][0]
   y0 = atoms.positions[i][1]
   z0 = atoms.positions[i][2]
   ro_fs_sum = 0
   ro_sc_sum = 0
   ro_gu_sum = 0
## Calculate energy of pair atoms
   for l in range(0, h, 1):
       rij = sqrt((x0-temp_clust[l][0])**2 + (y0-temp_clust[l][1])**2 + (z0-temp_clust[l][2])**2)
       if rij != 0:
          ro_fs_i = A*A*((rij-d)**2 + B*B*(rij-d)**4)
          ro_fs_sum = ro_fs_sum + ro_fs_i
          ro_sc_i = (a/rij)**m
          ro_sc_sum = ro_sc_sum + ro_sc_i
          ro_gu_i = exp(-2*q*((rij-r0)/r0))
          ro_gu_sum = ro_gu_sum + ro_gu_i
       ro_fs.append(ro_fs_sum)
       ro_sc.append(ro_sc_sum)
       ro_gu.append(ro_gu_sum)
## write to file
file1 = open('Surface properties_1K.xyz', 'w')
file1.write('%s\n\n' % surf_nn)
for i in range(0, nn, 1):
   if int(tags[i]) == 1:
      file1.write('%s\t' % atoms.symbols[i])
      file1.write('%s\t'  % atoms.positions[i][0])
      file1.write('%s\t'  % atoms.positions[i][1])
      file1.write('%s\t'  % atoms.positions[i][2])
      file1.write('%s\t' % surf_energy[i])
      file1.write('%s\t' % surf_energy_weight[i])
      file1.write('%s\t' % ro_fs[i])
      file1.write('%s\t' % ro_sc[i])
      file1.write('%s\n' % ro_gu[i])
file1.close()


    
