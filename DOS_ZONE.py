from math import degrees, atan, sqrt

f = open('dos_opj.txt', 'r')
lines = f.readlines()
e = []
s = []
p = []
d = []
spd = []
while len(lines) > 0:
   line = lines.pop(0)
   e0, s0, p0, d0, spd0 = line.split()[:5]
   e.append(float(e0))
   s.append(float(s0))
   p.append(float(p0))
   d.append(float(d0))
   spd.append(float(spd0))
n = len(e)
s_pdos = 0
p_pdos = 0
d_pdos = 0
spd_pdos = 0
es_pdos = 0
ep_pdos = 0
ed_pdos = 0
espd_pdos = 0
ees_pdos = 0
eep_pdos = 0
eed_pdos = 0
eespd_pdos = 0
fs_pdos = 0
fp_pdos = 0
fd_pdos = 0
fspd_pdos = 0
for i in range(0, n, 1):
   s_pdos = s_pdos + s[i]
   p_pdos = p_pdos + p[i]
   d_pdos = d_pdos + d[i]
   spd_pdos = spd_pdos + spd[i]
   es_pdos = es_pdos + e[i] * s[i]
   ep_pdos = ep_pdos + e[i] * p[i]
   ed_pdos = ed_pdos + e[i] * d[i]
   espd_pdos = espd_pdos + e[i] * spd[i]
   ees_pdos = ees_pdos + e[i] * e[i] * s[i]
   eep_pdos = eep_pdos + e[i] * e[i] * p[i]
   eed_pdos = eed_pdos + e[i] * e[i] * d[i]
   eespd_pdos = eespd_pdos + e[i] * e[i] * spd[i]
   if e[i] <= 0:
      fs_pdos = fs_pdos + s[i]
      fp_pdos = fp_pdos + p[i]
      fd_pdos = fd_pdos + d[i]
      fspd_pdos = fspd_pdos + spd[i]
s_center = es_pdos / s_pdos
p_center = ep_pdos / p_pdos
d_center = ed_pdos / d_pdos
spd_center = espd_pdos / spd_pdos
s_width = sqrt(ees_pdos / s_pdos)
p_width = sqrt(eep_pdos / p_pdos)
d_width = sqrt(eed_pdos / d_pdos)
spd_width = sqrt(eespd_pdos / spd_pdos)
s_occ = fs_pdos / s_pdos
p_occ = fp_pdos / p_pdos
d_occ = fd_pdos / d_pdos
spd_occ = fspd_pdos / spd_pdos
file = open('dos_data.txt', 'w')
file.write('Orbitals\t' + 'zone center\t' + 'zone width\t' + 'zone occupation\n')
file.write('s orbital\t')
file.write('%s\t' % s_center)
file.write('%s\t' % s_width)
file.write('%s\n' % s_occ)
file.write('p orbital\t')
file.write('%s\t' % p_center)
file.write('%s\t' % p_width)
file.write('%s\n' % p_occ)
file.write('d orbital\t')
file.write('%s\t' % d_center)
file.write('%s\t' % d_width)
file.write('%s\n' % d_occ)
file.write('s+p+d orbital\t')
file.write('%s\t' % spd_center)
file.write('%s\t' % spd_width)
file.write('%s\n' % spd_occ)
   
