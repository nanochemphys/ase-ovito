from ase.io import read, write
from ase import Atom, Atoms

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
with open('cluster.xyz', 'w') as file:
    file.write(str(natoms) + '\n\n')
    for i in range(natoms):
        file.write(atom_symb[i] + '\t' + atom_x[i] + '\t' + atom_y[i] + '\t' + atom_z[i] + '\n')
file.close()
    
