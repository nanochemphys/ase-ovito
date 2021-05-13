from ase.io import read, write
from ase import Atom, Atoms
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier
from ovito.modifiers import ConstructSurfaceModifier

## Create surface layer by OVITO-tools
pipeline = import_file('cluster.xyz')
pipeline.modifiers.append(ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 3.1,
    select_surface_particles = True))
export_file(pipeline, 'clusters_tags.xyz', 'xyz', columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "Selection"])

##Read xyz-file with ovito-tags 
atoms = read('clusters_tags.xyz')
at_symb = atoms.symbols[0]
file1 = open('clusters_tags.xyz', 'r')
lines = file1.readlines()
while len(lines) > 0:
	tags = []
	symbols = []
	xx = []
	yy = []
	zz = []
	surf_nn = 0
	natoms = int(lines.pop(0))
	lines.pop(0)
	for _ in range(natoms):
		line = lines.pop(0)
		symbol, x, y, z, tag = line.split()[:5]
		surf_nn = surf_nn + int(tag)
		symbols.append(symbol)
		xx.append(float(x))
		yy.append(float(y))
		zz.append(float(z))
		tags.append(int(tag))
file1.close()
## write to out file
file2 = open('cluster_surf.xyz', 'w')
file2.write('%s\n\n' % surf_nn)
for i in range(0, natoms, 1):
	if int(tags[i]) == 1:
		file2.write('%s\t' % symbols[i])
		file2.write('%s\t'  % xx[i])
		file2.write('%s\t'  % yy[i])
		file2.write('%s\n'  % zz[i])
file2.close()


    
