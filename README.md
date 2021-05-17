# ase-ovito
scripts for calculate energy and electronic density for atomic cluster
Requirement programms
ASE 
https://wiki.fysik.dtu.dk/ase/install.html
OVITO
https://www.ovito.org/
NUMPY
https://numpy.org/install/

All this scripts are designed to get clusters data. Atomic structure of clusters specify using a xyz-file. Cluster data include surface energy and electron density in Finnis-Sinclair, Satton-Chen and Gupta potentials and MEAM model:

surface_energy.py
This script read input file - cluster_tmp.xyz, atomic structure of cluster, and make output file - xyz-file of surface atoms with additional columns. 
Two additional columns contain values of surface energy and surface energy weighted by the maximum coordination number.

semi.py
This script read input file - cluster_tmp.xyz, atomic structure of cluster, and make output file - xyz-file of surface atoms with additional columns. 
Three additional columns contain values of electron density, calculated in Finnis-Sinclair, Sutton-Chen and Gupta potentials.
Finnis-Sinclair potential:
https://doi.org/10.1088/0953-8984/18/19/008
Sutton-Chen potential:
https://doi.org/10.1016/j.jpcs.2015.03.008
Gupta potential:
https://doi.org/10.1080/10420150108211842

meam.py
This script read input file - cluster_tmp.xyz, atomic structure of cluster, and make output file - xyz-file of surface atoms with additional columns. 
Three additional columns contain values of electron density, calculated in MEAM-model by three different formulas
MEAM model:
doi:10.1016/j.susc.2003.12.043

Other scripts can be useful when working with atomic clusters:

neighborlist_ase.py
This script create list of neighbors for every cluster atom by ASE-tools with additional deleting probably dublicate atoms. Input file - xyz-file, output file - xyz-file with additional column with coordination number of the corresponding atom

surface_ovito.py
This script create xyz-files: one with additional column of tags, 1 - suface atom, 0 - volume atom, another with only surface atoms.

DOS_ZONE.py
This script calculate centre, widths and occupation s, p, d-orbitals density of state

revcon_xyz.py
REVCON is DL_POLY file, contained molecule dynamic information. This script transform REVCON to xyz-file
