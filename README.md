The above scripts are used to calculate the surface energy in the EMT potential, the electron density in the Gupta, Sutton-Chen, Finnis-Sinclair, and MEAM potentials.
The scripts are launched with the command:
python3 name_of_script.py
The input file contains the coordinates of the cluster in xyz format or REVCON (output file DL_Poly).
The following libraries are required for the scripts to work:
Ovito https://www.ovito.org/
ASE https://wiki.fysik.dtu.dk/ase/
NUMPY https://numpy.org/
Development and testing of codes was carried out on operating systems CentOS 7 and Ubuntu …. 

Authors: Nadezhda V. Dokhlikova (dohlikovanv@gmail.com), Sergey V. Doronin (sedoronin@gmail.com)
Request when using scripts to quote the publication:
Sergey V. Doronin, Nadezhda V. Dokhlikova, Maxim V. Grishin Descriptor of catalytic activity nanoparticles surface: Atomic and molecular hydrogen on gold // …

1. Description of scripts
1.1. surface_energy.py
The calculation of the surface energy of the selected atom in the nearest environment is carried out in the potential EMT
The output file name_of_file.xyz is a list of surface atoms with their symbols and coordinates given in the first four columns. The fifth shows the number of atoms from the nearest environment of the surface atom. The sixth and seventh columns show the surface energy and its value weighted by the maximum coordination number for each corresponding surface atom.
The general algorithm can be broken down into several main blocks:
1) The search for surface atoms is carried out using the tools of the OVITO program, the cutoff radius is 3.1 A.
2) For each found surface atom the nearest environment is constructed using ASE tools, the radius of determining the nearest environment is 2.3 * covalent radius.
3) In the constructed environment the surface energy of the selected atom is calculated
E_i=1/N*∑(E_ij) , where E_ij – the energy of a pair of i- and j- atoms, selected from the immediate environment. After summing up all pairs over the environment, the energy of the i-th atom is normalized to the number of atoms in the environment (N).
1.2. semi.py
The script calculates the electron density of the selected atom in the immediate environment in the Gupta, Sutton-Chen, Finnis-Sinclair potentials.
The output file name_of_output.xyz is a list of surface atoms with their symbols and Cartesian coordinates given in the first four columns. In the fifth, the number of atoms from the nearest environment of the surface atom is given. The sixth, seventh and eighth columns show the values of the electron density of atoms in the Finnis-Sinclair, Sutton-Chen and Gupta potentials, respectively.
The general algorithm consists of the following blocks:
1. The search for surface atoms is carried out using the tools of the OVITO program, the cutoff radius is 3.1 A
2. For each found surface atom the nearest environment is constructed using the ASE tools, the radius of determining the nearest environment is 2.3 * the covalent radius of the atom.
3. Taking into account the nearest environment, electron densities are calculated for each surface atom.
The electron densities at the indicated potentials are calculated as follows:
- Gupta potential: https://doi.org/10.1080/10420150108211842
- Sutton-Chen potential: https://doi.org/10.1016/j.jpcs.2015.03.008
- Finnis-Sinclair potential: https://doi.org/10.1088/0953-8984/18/19/008
1.3. meam.py
Meam.py calculates the electron density of the selected atom in the MEAM model.
The output file name_of_output.xyz is a list of surface atoms, in which their symbols and Cartesian coordinates are given in the first four columns. In the fifth column, the number of the nearest neighbor atoms of the surface atom is given. The sixth, seventh and eighth columns show the values of the electron density of atoms in the MEAM model in various approximations.
The general script algorithm consists of the following blocks:
1. The search for surface atoms is carried out using the tools of the OVITO program ({method Alpha-shape}) with probe sphere radius 3.1 Å.
2. For each surface atom the neighbor list is constructed using the ASE tools, the radius of determining the neighbor atoms is 2.3 * the covalent radius of the atom.
3. Taking into account the neighbor list, electron densities are calculated for each surface atom.
The MEAM model was described in the text of the paper doi:10.1016/j.susc.2003.12.043.

2. Supporting scripts
2.1. revcon_xyz.py
Script converts the REVCON file created in the DL_POLY molecular modeling software package to the xyz format. The output file is written as cluster.xyz
2.2. surface_ovito.py
Surface_ovito.py is used to search for surface atoms of cluster using the Alpha-shape method with probe sphere radius equal to 4 and smoothing level equal to 8 integrated into the OVITO program. The script processes the file name_of_input.xyz with the atomic structure of the cluster and creates two output files clus-ter_tags.xyz with all the atoms of the cluster being processed, containing an additional column of tags, where 1 is an atom on the cluster surface, 0 is an atom in the volume of the cluster, and cluster_surf.xyz containing only surface atoms.
2.3. dos_zone.py
The script calculates the characteristics of the density of states of atoms, center, width and population based on DOS data.
The input file is a DOS file calculated in Quantum Espresso, containing 5 columns: the level of electronic energy (must be corrected for the Fermi level), the densities of states of s-, p-, d-bands and their sum. The output file contains four columns: band type, center, width and band filling.
2.4. neighborlist_ase.py
The script builds neighbor list of atoms using ASE tools with the removal of possible duplicates. During a check of the ASE library, it was found that the Neighbourlist function in some cases re-writes some atoms in the neighbor list. This point has been corrected in the current script.
The input file is the xyz file. The output file, in addition to the standard xyz output, contains an additional column with a coordination number for each atom. Such a script can be useful when constructing coordination numbers for atoms or when checking various scripts in which various parameters are calculated based on the neighbor list.
