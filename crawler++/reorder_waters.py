import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

system = mda.Universe("Li+W20.xyz")
watersO = system.select_atoms("name O")
watersH = system.select_atoms("name H")
pos = watersH.positions.copy()
# for atom in system.atoms:
#     print(atom.position)
# exit()
r1 = system.atoms[1].position - system.atoms[0].position
r2 = system.atoms[5].position - system.atoms[0].position
# r1 /= np.linalg.norm(r1)
# r2 /= np.linalg.norm(r2)
print(np.linalg.norm(r1), np.linalg.norm(r2))

# print((180 / np.pi) * np.arccos(np.dot(r1, r2)))

oxygensDict = {}
for i in range(len(watersO)):
    dists = distance_array(watersO[i].position, pos)[0]
    print(i)
    print(dists)
    oxygensDict[i] = []
    for j in range(len(dists)):
        # print(j)
        if dists[j] < 1.01:
            oxygensDict[i].append(watersH[j].ix)
# check = []
# for _, value in oxygensDict.items():
#     check += [value[0], value[1]]
# print(len(check))
# print(len(set(check)))
# exit()
newAtomGroup = []
for i in range(len(watersO)):
    newAtomGroup += [system.atoms[watersO[i].ix], system.atoms[oxygensDict[i][0]], system.atoms[oxygensDict[i][1]]]

# newAtomGroup += [system.atoms[-1]]
newAtomGroup = mda.AtomGroup(newAtomGroup)
newAtomGroup.write("Li_waterClusterReordered.xyz")

system = mda.Universe("Li_waterClusterReordered.xyz")
# for atom in system.atoms:
#     print(atom)
r1 = system.atoms[4].position - system.atoms[3].position
r2 = system.atoms[5].position - system.atoms[3].position
print(np.linalg.norm(r1), np.linalg.norm(r2))
exit()
r1 /= np.linalg.norm(r1)
r2 /= np.linalg.norm(r2)

print((180 / np.pi) * np.arccos(np.dot(r1, r2)))
