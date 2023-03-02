from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
import pandas as pd
from sys import argv
init()

pose = pose_from_pdb(argv[1])

total_resi = pose.total_residue()

sf = get_fa_scorefxn()
sf(pose)

df = pd.DataFrame(columns=['residueNo','resi_En','chain','fa_elec'])

resils = []
Enls = []
chainls = []
elecls = []
for i in range(1, total_resi+1):
	resi_info = pose.pdb_info().pose2pdb(i)
	resi_info = resi_info.rsplit(' ')
	resi = resi_info[0]
	chain = resi_info[1]
	
	resils.append(resi)
	chainls.append(chain)
	En = pose.energies().residue_total_energy(i)
	Enls.append(En)

	elec = pose.energies().residue_total_energies(i)[fa_elec]
	elecls.append(elec)

df.residueNo = resils
df.resi_En = Enls
df.chain = chainls
df.fa_elec = elecls

df.to_csv(argv[1][:-4]+'.csv',index=False)



