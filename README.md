# ase-ovito
scripts for calculate energy and electronic density for atomic cluster
Requirement programms
ASE 
https://wiki.fysik.dtu.dk/ase/install.html
OVITO
https://www.ovito.org/
NUMPY
https://numpy.org/install/

surface_energy.py
semi.py
meam.py
This scripts read input file - cluster_tmp.xyz, atomic structure of cluster, and make output file - xyz-file of surface atoms with additional columns. 
In surface_energy.py columns contain values of surface energy and surface energy weighted by the maximum coordination number.
In semi.py columns contain values of electron density, calculated in Finnis-Sinclair, Satton-Chen and Gupta potentials.
In meam.py columns contain values of electron density, calculated in MEAM-model

neighborlist_ase.py
This script create list of neighbors for every cluster atom by ASE-tools with additional deleting probably dublicate atoms. Input file - xyz-file, output file - xyz-file with additional column with coordination number of the corresponding atom

surface_ovito.py
This script create xyz-files: one with additional column of tags, 1 - suface atom, 0 - volume atom, another with only surface atoms.

DOS_ZONE.py
This script calculate centre, widths and occupation s, p, d-orbitals density of state

revcon_xyz.py
REVCON is DL_POLY file, contained molecule dynamic information. This script transform REVCON to xyz-file
