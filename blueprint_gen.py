from sys import argv
from pyrosetta import *
from pyrosetta.toolbox import *
init()

input_stuct = argv[1]

pose = pose_from_pdb(input_stuct)
sq = pose.sequence()


f = open(input_stuct[:-4]+'.blueprint', 'a+')
for i in range(len(sq)):
	f.write(str(i+1)+' '+sq[i]+' . \n')
f.close()