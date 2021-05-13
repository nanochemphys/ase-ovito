from ase.io import read, write
from ase import Atom, Atoms
from ase.calculators.emt import EMT
from ase import neighborlist
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.data import covalent_radii
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier
from ovito.modifiers import ConstructSurfaceModifier
import numpy as np
from math import sqrt, exp

###############################
#########SUBROUTINGS###########
###############################

##SUBRUTINE r_ij
def r_i_j(r_i, r_j): return sqrt((r_j[0] - r_i[0])**2 + (r_j[1] - r_i[1])**2 + (r_j[2] - r_i[2])**2)

## SUBRUTINE Vector Commponent
def proj(alpha, r_i, r_j): 
    if alpha == 0:
        return sqrt((r_j[0] - r_i[0])**2)/r_i_j(r_i, r_j)
    elif alpha == 1:
        return sqrt((r_j[1] - r_i[1])**2)/r_i_j(r_i, r_j)
    else:
        return sqrt((r_j[2] - r_i[2])**2)/r_i_j(r_i, r_j)
        
##Partial electron density 3-ORDER
def ro_i_3(r_i, r_j):
    summ1 = 0
    for a in 0, 1, 2:
        for b in 0, 1, 2:
            for g in 0, 1, 2:
                for j, r_jj in enumerate(r_j):
                    r_jj = r_j.pop(j)
                    if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
                        tmp = r_jj
                    else:
                        summ1 = summ1 + proj(a, r_i, r_jj) * proj(b, r_i, r_jj) * proj(g, r_i, r_jj) * exp(-1* b_3*(r_i_j(r_i, r_jj)/r_e -1))
                        r_j.insert(j, r_jj)
    summ2 = 0
    for a in 0, 1, 2:
        for j, r_jj in enumerate(r_j):
            r_jj = r_j.pop(j)
            if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
                tmp = r_jj
            else:
                summ2 = summ2 + proj(a, r_i, r_jj) * exp(-1* b_3*(r_i_j(r_i, r_jj)/r_e -1))
                r_j.insert(j, r_jj)
    return summ1 - 3 * summ2 / 5
    
## Partial electron density 2-ORDER
def ro_i_2(r_i, r_j):
    summ1 = 0
    for a in 0, 1, 2:
        for b in 0, 1, 2:
            for j, r_jj in enumerate(r_j):
                r_jj = r_j.pop(j)
                if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
                    tmp = r_jj
                else:
                    summ1 = summ1 + proj(a, r_i, r_jj) * proj(b, r_i, r_jj) * exp(-1* b_2*(r_i_j(r_i, r_jj)/r_e -1))
                    r_j.insert(j, r_jj)
    summ2 = 0
    for j, r_jj in enumerate(r_j):
        r_jj = r_j.pop(j)
        if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
            tmp = r_jj
        else:
            summ2 = summ2 + exp(-1* b_2*(r_i_j(r_i, r_jj)/r_e -1))
            r_j.insert(j, r_jj)
    return summ1 - summ2/3
    
## Partial electron density 1-ORDER
def ro_i_1(r_i, r_j):
    summ = 0
    for a in 0, 1, 2:
        for j, r_jj in enumerate(r_j):
            r_jj = r_j.pop(j)
            if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
              tmp = r_jj
            else:
                summ = summ + proj(a, r_i, r_jj) * exp(-1* b_1*(r_i_j(r_i, r_jj)/r_e -1))
                r_j.insert(j, r_jj)
    return summ

## Partial electron density 0-ORDER
def ro_i_0(r_i, r_j):
   summ = 0
#   print('i-atom')
#   print(r_i)
   for j, r_jj in enumerate(r_j):
      r_jj = r_j.pop(j)
      if r_i[0] == r_jj[0] and r_i[1] == r_jj[1] and r_i[2] == r_jj[2]:
          tmp = r_jj
      else:
          summ = summ + exp(-1* b_0*(r_i_j(r_i, r_jj)/r_e -1))
#      print(screen(r_i, r_jj, r_j))
          r_j.insert(j, r_jj)
#      print(summ)
   return summ

## Gamma function
def gamma_i(r_i, r_j): return t_1 * ((ro_i_1(r_i, r_j)/ro_i_0(r_i, r_j)))**2 + t_2 * ((ro_i_2(r_i, r_j)/ro_i_0(r_i, r_j)))**2 + t_3 * ((ro_i_3(r_i, r_j)/ro_i_0(r_i, r_j)))**2

## electron density
def el_dens_meam(r_i, r_j):
   ro_i_g1 = ro_i_0(r_i, r_j) * sqrt(1 + gamma_i(r_i, r_j)) ## I
   try:
       ro_i_g2 = ro_i_0(r_i, r_j) * exp(gamma_i(r_i, r_j) / 2) ## II
   except OverflowError:
       ro_i_g2 = float('inf')
   ro_i_g3 = ro_i_0(r_i, r_j) * 2 / (1+exp(-1 * gamma_i(r_i, r_j))) ## III
   return ro_i_g1, ro_i_g2, ro_i_g3

########################################
###############MAIN#####################
########################################

##############ovito##############################
pipeline = import_file('cluster_tmp.xyz')
pipeline.modifiers.append(ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 3.1,
    select_surface_particles = True))
# Export results of the clustering algorithm to a text file:
export_file(pipeline, 'clusters_out.xyz', 'xyz', columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "Selection"])

#############Read cluster xyz file####################
atoms = read('clusters_out.xyz')
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
nn = atoms.get_global_number_of_atoms()
if atoms.symbols[0] == 'Cu':
	t_1 = 3.14
    t_2 = 2.49
    t_3 = 2.95
    b_0 = 3.63
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.56
elif  atoms.symbols[0] == 'Ag':
    t_1 = 5.54
    t_2 = 2.45
    t_3 = 1.29
    b_0 = 4.46
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.88
elif atoms.symbols[0] == 'Au':
    t_1 = 1.59
    t_2 = 1.51
    t_3 = 2.61
    b_0 = 5.45
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.88
elif atoms.symbols[0] == 'Ni':
    t_1 = 3.57
    t_2 = 1.60
    t_3 = 3.70
    b_0 = 2.45
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.49
elif atoms.symbols[0] == 'Pd':
    t_1 = 2.34
    t_2 = 1.38
    t_3 = 4.48
    b_0 = 4.98
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.75
elif atoms.symbols[0] == 'Pt':
    t_1 = 2.73
    t_2 = -1.38
    t_3 = 3.29
    b_0 = 4.67
    b_1 = 2.2
    b_2 = 6.0
    b_3 = 2.2
    r_e = 2.77
## create ase-neighborlist, natural cuttoffs are multiple covalent radii
cutoff = neighborlist.natural_cutoffs(atoms, 2.3)
nl = NeighborList(cutoff, skin=0, self_interaction=True, bothways=True)
nl.update(atoms)
at_symb = atoms.symbols[0]
with open('%s_clust%s_surf%s.xyz' % (at_symb, nn, surf_nn), 'w') as file:
    file.write('%s\n' % surf_nn)
    file.write('Symbol\t' + 'x_I\t\t' + 'y_I\t\t' + 'z_I\t\t' + 'number_I\t\t' + 'ro_I\t\t' + 'ro_II\t\t' + 'ro_III\t\t' + 'ro_FinSin\t\t' + 'ro_SutChen\t\t' + 'ro_Gupta\n')
    for i, atom in enumerate(atoms):
       if tags[i] == 1:
           indices, offsets = nl.get_neighbors(i)
           temp_clust = []
           for j, offset in zip(indices, offsets):
                temp_clust.append(atoms.positions[j]) ## append to list
           atoms_i = atoms.positions[i]
           temp_clust = [temp for temp in temp_clust if (temp[0] != atoms_i[0] and temp[1] != atoms_i[1] and temp[2] != atoms_i[2])]
           len_nn = len(temp_clust)
           ro_i_g1, ro_i_g2, ro_i_g3 = el_dens_meam(atoms_i, temp_clust)
           file.write('%s\t' % atoms.symbols[i])
           file.write('%s\t' % atoms.positions[i][0])
           file.write('%s\t' % atoms.positions[i][1])
           file.write('%s\t' % atoms.positions[i][2])
           file.write('%s\t' % len_nn)
           file.write('%s\t'  % ro_i_g1)
           file.write('%s\t'  % ro_i_g2)
           file.write('%s\n'  % ro_i_g3)
file.close()
