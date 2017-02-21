import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style("whitegrid")

from msmbuilder import dataset

# load trajectories

Abl_in = dataset.MDTrajDataset("../Src_and_Abl_pdbs/Abl_DFG_in/ABL1*.pdb")
Abl_out = dataset.MDTrajDataset("../Src_and_Abl_pdbs/Abl_DFG_out/ABL1*.pdb")

Src_in = dataset.MDTrajDataset("../Src_and_Abl_pdbs/Src_DFG_in/SRC*.pdb")
Src_out = dataset.MDTrajDataset("../Src_and_Abl_pdbs/Src_DFG_out/SRC*.pdb")

# Here we are defining anew a coordinate to define the DFG flip.
# Okay this is my last attempt before going into RMSD or scikit learn.
# In this case we are trying the distance between the
# CZ of the F in the DFG motif and the CA of a stable glycine/alanine.
# QIAS_G_MAY in Src and QISS_A_MEY in Abl
# Below are these coordinates in PDBs without hydrogens
# These are atom numbers acquired crudely by looking at the PDB.
Abl_DFG_distance = [871,1136]
Src_DFG_distance = [835,1093]

def DFG_distance(trajectories,def_DFG):

    distance = []

    for traj in trajectories:

        topology = traj.topology
        print 'Atom distances computed between %s and %s' %(topology.atom(def_DFG[0]),topology.atom(def_DFG[1]))        

        distance.append(md.compute_distances(traj,[def_DFG]))

    flattened_distance = np.asarray([val for sublist in distance for val in sublist])

    return [flattened_distance]

[Abl_in_distance] = DFG_distance(Abl_in, Abl_DFG_distance)
[Abl_out_distance] = DFG_distance(Abl_out, Abl_DFG_distance)
[Src_in_distance] = DFG_distance(Src_in, Src_DFG_distance)
[Src_out_distance] = DFG_distance(Src_out, Src_DFG_distance)

plt.figure(figsize=(10,3))

for i in range(len(Abl_in_distance)):
    if i == 0:
        plt.scatter(Abl_in_distance[i],0.5, edgecolors="r", marker='o', linewidth='3', s=80, facecolors='none',label='ABL DFG-in')
    else:
        plt.scatter(Abl_in_distance[i],0.5, edgecolors="r", marker='o', linewidth='3', s=80, facecolors='none')
for i in range(len(Abl_out_distance)):
    if i == 0:
        plt.scatter(Abl_out_distance[i],0.25, edgecolors="m", marker='o', linewidth='3', s=80, facecolors='none',label='ABL DFG-out')
    else:
        plt.scatter(Abl_out_distance[i],0.25, edgecolors="m", marker='o', linewidth='3', s=80, facecolors='none')
for i in range(len(Src_in_distance)):
    if i == 0:
	plt.scatter(Src_in_distance[i],-0.25, edgecolors="b", marker='o', linewidth='3', s=80, facecolors='none',label='SRC DFG-in')
    else:
        plt.scatter(Src_in_distance[i],-0.25, edgecolors="b", marker='o', linewidth='3', s=80, facecolors='none')
for i in range(len(Src_out_distance)):
    if i == 0:
        plt.scatter(Src_out_distance[i],-0.5, edgecolors="g", marker='o', linewidth='3', s=80, facecolors='none',label='SRC DFG-out')
    else:
        plt.scatter(Src_out_distance[i],-0.5, edgecolors="g", marker='o', linewidth='3', s=80, facecolors='none')


plt.xlabel('Distance (nm)')
plt.ylim((-0.75, 0.75))
plt.yticks([])
plt.legend()

plt.tight_layout()
plt.savefig('DFG_distance_1D.png',dpi=300)


