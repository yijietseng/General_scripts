from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
import sys,glob
from sys import argv
init()

if '-h' in argv:
	Description = '''
	tag			argument		Description
	-m			<No_arg>		Execute multi-structure measurement. This method can only do structure names in series 
	-l			<No_arg>		Measure the SC of the ligand against the complex
	-r			<path>			To check if record the SC into a log. Argument is the path of the output logfile
	-p			<path>			input the PDB structure. For multiple structure mode, input only the path to the structure directory
	'''
	print(Description)
	sys.exit()


if '-p' in argv:
	PDBname = argv[argv.index('-p')+1]

if '-m' in argv:
	if '-r' in argv:
		if '-p' not in argv:
			print('Please input a path to the PDB structure!!!')
		else:
			strucNo = argv.index('-m')+1
			logpath = argv[argv.index('-r')+1]
			PDB = glob.glob(PDBname+'*.pdb')
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\n')
			sc = ShapeComplementarityFilter()
			SCls = []
			for i in PDB:
				pose = pose_from_pdb(i)
				sc_score = sc.score(pose)
				f.write(i+'\t'+str(sc_score)+'\n')
				SCls.append(i+':'+str(sc_score))
			f.close()
			print(SCls)
	else:
		if '-p' not in argv:
			print('Please input a path to the PDB structure!!!')
		else:
			strucNo = argv.index('-m')+1
			logpath = argv[argv.index('-r')+1]
			PDB = glob.glob(PDBname+'*.pdb')
			sc = ShapeComplementarityFilter()
			SCls = []
			for i in PDB:
				pose = pose_from_pdb(i)
				sc_score = sc.score(pose)
				SCls.append(i+':'+str(sc_score))
			print(SCls)
else:
	if '-p' not in argv:
		print('Please input a path to the PDB structure!!!')
	else:
		sc = ShapeComplementarityFilter()
		pose = pose_from_pdb(PDBname)
		sc_score = sc.score(pose)
		print('SC score for '+PDBname+' is:'+str(sc_score))
if '-l' in argv:
	if '-p' not in argv:
		print('Please input a path to the PDB structure!!!')
	else:
		pose = pose_from_pdb(PDBname)
		
		total = pose.total_residue()

		sc = ShapeComplementarityFilter()
		
		ligstr=''
		enzstr=''
		for i in range(1,total+1):
			resi_info = pose.pdb_info().pose2pdb(i)
			resi = resi_info.replace(' ','')
			if resi[-1] == 'A':
				ligstr += resi+','
			else:
				enzstr += resi+','
		ligstr = ligstr[:-1]
		enzstr = enzstr[:-1]
		
		sc.residues1(ligstr)
		sc.residues2(enzstr)
		sc_score = sc.score(pose)
		print('SC score for '+PDBname+' is:'+str(sc_score))


