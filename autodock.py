'''
Created by Yi Jie (Josh) Tseng --- 2/24/2022
Modified by Yi Jie (Josh) Tseng --- 9/2/2023

Software requirement: PyMOL and Python3

Description:
This script is used to prepare receptor and ligands for input, define search space and run autodock vina. The output of results and analysis is still under development. Before running autodock, few things are required: 
	1. A search box and a center position
	2. PDBQT file of reveptor 
	3. PDBQT file of ligand(s)
The instruction below will go over how to prepare for those requirements. 

Instructions:

1. Set up working directory:
	This script can create folders to organize inputs, outputs and logs. To do this, in a terminal, type python3 autodock.py -s <XXX> True -f <No.>. -s flag is to set up folders. The first argument is the folder name, the "True" in the second argument is to turn on or off the numbering mechanism. -f option is to specify a certain working number. Please see below examples:
	1.1 Specify folder name + numbering mechanism on. If this command is executed multiple times, it will generate folders with specified name and ordered two digit numbers. But the numbering mechanism can be overwritten by the -f flag (see 4th example).
		Ex:
			First: python3 autodock.py -s test True		->	test01/  
			Second: python3 autodock.py -s test True	->	test02/
			Third:  python3 autodock.py -s test True	->	test03/
			Fourth:  python3 autodock.py -s test True -f 8	->	test08/

	1.2 Specify folder name + numbering mechanism off. -f flag is recommended. If the command is executed multiple times without -f, it could potentially overwrite a previous work. 
		Ex:
			First: python3 autodock.py -s test	->	test/  
			Second: python3 autodock.py -s test	->	test/  (Overwritting the first run)
			Third: python3 autodock.py -s test -f 1  ->   test01

	1.3 If no name specified, by default, the folder name is "Dock" and the numbering mechanism is automatically on. Again the numbering mechanism can be overwritten by -f flag.
		Ex:
			First: python3 autodock.py -s   ->   Dock01/
			Second: python3 autodock.py -s	->	Dock02/
			Third:  python3 autodock.py -s -f 5	->	Dock05/


2. Identifying seach box and center position:
	****IMPORTANT****: PyMOL is required to visualize the search box and center position in this script. 
	
	If the binding site of the target is known, select a set of binding site residues. This script is able to calculate the center coordinace of the selected residues based on the average of the CA atom coordinances of the selected residues. If only a single residue is selected, for example, a bound ligand. By selecting the ligand molecule, it can calculate the cooridinace of the center position of the molecule based on the average of the cooridinance of all atoms. If one single atom is selected, it will use the coordinance of the selected atom as the center position. To do the above steps, navegate to the working directory in PyMOL, run this script in PyMOL by typing <run autodock.py>, select residue(s) by clicking on the desired residues (show sticks representation might be required) and in the command line box type Box('sele', X, Y, Z), where 'sele' is the seleted residues; X, Y and Z are the dimention of the box. Once the center position and the dimention of the search box is optimal, save the configuration file running the same function by typing Box('sele', X, Y, Z, save='WORKING_FOLDER/config.txt'), the name of the configuration can be arbitrary. But it is recommended to include the name of the working folder.

	
3. Preparation of receoptor and ligand PDBQT files:
	To prepare the required PDBQT files for receptor and ligands (XXXX.pdbqt), through passing -p flag, this script runs scripts from mgltools, which is a software for molecular structure analysis and can convert PDB files to PDBQT files and assign partial charges to all atoms. if the initial working folder set up does not run with this process at the same time, a -w flag specifying the working directory is required. To do this, run the script in the terminal by typing:
		With initial working directory set up:
			python3 autodock.py -s -p 
		Without initial working directory set up:
			python3 autodock.py -p -w Dock01 -c Dock01/config.txt 

			
4. Run autodock:
	Before running autodock, receptorPDBQT, ligandPDBQT and configuration file are required. Like in step 3, if the initial working directory set up is not executed at the same time, it will also require -w for specifying working directory. To run docking, -d flag is passed. 
		Ex:
			python3 autodock.py -d -r Dock01/Receptor/receptor.pdbqt -l Dock01/Ligand/ligand.pdbqt -c Dock01/config.txt -w Dock01


7. Analyze data:
	Under development


'''

import os, itertools, argparse, glob, sys, subprocess



def center_coord(selection:str):
	import pymol
	from pymol import cmd
	
	#--------------------------------------IMPORTANT--------------------------------------
	# This function can only be run in pymol

	# get selection residue list
	pymol.stored.resi_ls = []
	cmd.iterate(f'({selection})', 'stored.resi_ls.append(resi)')
	# remove duplicates and sort
	pymol.stored.resi_ls = list(set(pymol.stored.resi_ls))

	# get chain ID of the selected residues
	chainID = cmd.get_chains(selection)[0]

	# calculate the center coordinance. If a group of residues are selected, it will calculate the averaged coordinance based on the CA of the selected residues.  
	if len(pymol.stored.resi_ls) > 1:
		# get all CA cooridinances from the selection
		ca_coords = []
		#print(pymol.stored.resi_ls)
		for resi in pymol.stored.resi_ls:
			#print(f'chain {chainID} and resi {resi} and name CA')
			cmd.select(f'ca_atoms', f'chain {chainID} and resi {resi} and name CA')
			if cmd.get_coords('ca_atoms',1) is not None:
				ca_coords.append(list(cmd.get_coords('ca_atoms',1)[0]))

		#print(ca_coords)
		# calculate center cooridinance by averaging all selected CA
		total_x = total_y = total_z = 0
		for coord in ca_coords:
			x, y, z = coord
			total_x += x
			total_y += y
			total_z += z
		cen_coord = [total_x/len(ca_coords), total_y/len(ca_coords), total_z/len(ca_coords)]

		return cen_coord
	
	# If a single residue is selected, it will calculate the averaged coordinance based on all atoms.  
	elif len(pymol.stored.resi_ls) == 1:
		pymol.stored.atom_ls = []
		cmd.iterate(f'{selection}', 'stored.atom_ls.append(name)')
		
		# if only one atom is selected, it will use the coordinance of the selected atom
		if len(pymol.stored.atom_ls) == 1:
			cen_coord = list(cmd.get_coords(selection,1)[0])
			
			return cen_coord
		else:
			atom_coords = []
			for atom in pymol.stored.atom_ls:
				cmd.select('atom', f'{selection} and name {atom}')
				atom_coords.append(list(cmd.get_coords('atom',1)[0]))
			
			total_x = total_y = total_z = 0
			for coord in atom_coords:
				x, y, z = coord
				total_x += x
				total_y += y
				total_z += z
			cen_coord = [total_x/len(atom_coords), total_y/len(atom_coords), total_z/len(atom_coords)]

			return cen_coord

	else:
		raise ValueError('Invalid selection!!!')


def Box(selection:str, x:float, y:float, z:float, save=''):
	'''
	Run the script in pymol to visualize the search area.
	Sets up the search box within the protein, which is
	used in the docking protocol.
	'''
	#--------------------------------------IMPORTANT--------------------------------------
	# This function can only be run in pymol
	import pymol
	from pymol import cmd
	from pymol.cgo import LINEWIDTH, BEGIN, LINES, VERTEX, END

	
	# Clear box and center point from a previous run
	if 'Position' in cmd.get_object_list():
		cmd.delete('Position')
		cmd.delete('Box')

	# Get center coordinance
	pX, pY, pZ = center_coord(selection)
	
	# create a pseudoatom for visualization
	cmd.pseudoatom('Position', pos=[pX, pY, pZ])
	cmd.show('spheres', 'Position')
	#cmd.center('Position')
	
	# calculate the minimum and maximum of XYZ based on the identified size
	minX = pX-float(x)
	minY = pY-float(y)
	minZ = pZ-float(z)
	maxX = pX+float(x)
	maxY = pY+float(y)
	maxZ = pZ+float(z)

	'''
	   4-----------8
	  /|          /|
	 / |         / |
	3-----------7  |				↑
	|  |    O   |  |				|   
	|  2--------|--6				Y  /
	| /         | /					| Z
	|/          |/					|/
	1-----------5          ←----X---o
	
	This is a schematic of the vertices and lines of the box. The codes below follow the same numbering scheme. "O" represents the center of the box.
	'''

	boundingBox = [
		LINEWIDTH, float(6),
		BEGIN, LINES,
		VERTEX, minX, minY, minZ, VERTEX, minX, minY, maxZ, # drawing line between #1-2
		VERTEX, minX, maxY, minZ, VERTEX, minX, maxY, maxZ, # drawing line between #3-4
		VERTEX, maxX, minY, minZ, VERTEX, maxX, minY, maxZ, # drawing line between #5-6
		VERTEX, maxX, maxY, minZ, VERTEX, maxX, maxY, maxZ, # drawing line between #7-8

		VERTEX, minX, minY, minZ, VERTEX, maxX, minY, minZ, # drawing line between #1-5
		VERTEX, minX, maxY, minZ, VERTEX, maxX, maxY, minZ, # drawing line between #3-7
		VERTEX, minX, maxY, maxZ, VERTEX, maxX, maxY, maxZ, # drawing line between #4-8
		VERTEX, minX, minY, maxZ, VERTEX, maxX, minY, maxZ, # drawing line between #2-6

		VERTEX, minX, minY, minZ, VERTEX, minX, maxY, minZ, # drawing line between #1-3
		VERTEX, maxX, minY, minZ, VERTEX, maxX, maxY, minZ, # drawing line between #5-7
		VERTEX, minX, minY, maxZ, VERTEX, minX, maxY, maxZ, # drawing line between #2-4
		VERTEX, maxX, minY, maxZ, VERTEX, maxX, maxY, maxZ, # drawing line between #6-8
		END]
	boxName = 'Box'

	# load drawn box
	cmd.load_cgo(boundingBox, boxName)

	if save:
		with open(save, 'w') as fout:
			fout.write(f'center_x = {pX}\n')
			fout.write(f'center_y = {pY}\n')
			fout.write(f'center_z = {pZ}\n')
			fout.write(f'size_x = {x}\n')
			fout.write(f'size_y = {y}\n')
			fout.write(f'size_z = {z}\n')


def setup_folders(folder_name='', file_incre='', numbering=False):
	'''
	This function can help set up necessary folders. If folder_name is not given, by default, the folder name is set to DockNN/, where NN is the two digit number. A user can not only specify the folder name, but can also specify whether a numbering system is turn on. 
	'''
	if folder_name:
		if numbering:
			# Number the folders or specify a number
			if file_incre:
				file_incre = file_incre
			else:
				for i in itertools.count(start=1):
					if not os.path.exists(f'{folder_name}{i:02}/'):
						file_incre = f'{i:02}'
						break    		
			try:
				Folder_ls = [f'{folder_name}{file_incre}', f'{folder_name}{file_incre}/Receptors', f'{folder_name}{file_incre}/Ligands', 
							f'{folder_name}{file_incre}/output', f'{folder_name}{file_incre}/logs']
				for folder in Folder_ls:
					os.mkdir(folder)
				print('Folders created!!!1')
			except:
				print('Folder exist!!!1')
				pass

			return f'{folder_name}{file_incre}'
		else:
			if file_incre:
				try:
					Folder_ls = [f'{folder_name}{file_incre}', f'{folder_name}{file_incre}/Receptors', f'{folder_name}{file_incre}/Ligands', 
								f'{folder_name}{file_incre}/output', f'{folder_name}{file_incre}/logs']
					for folder in Folder_ls:
						os.mkdir(folder)
					print('Folders created!!!2')
				except:
					print('Folder exist!!!2')
					pass
			else:
				try:
					Folder_ls = [f'{folder_name}', f'{folder_name}/Receptors', f'{folder_name}/Ligands', 
								f'{folder_name}/output', f'{folder_name}/logs']
					for folder in Folder_ls:
						os.mkdir(folder)
					print('Folders created!!!3')
				except:
					print('Folder exist!!!3')
					pass

			return folder_name
	else:
		# if no folder_name given, the numbering system is automatically on
		numbering = True
		if file_incre:
			file_incre = file_incre[0]
		else:
			for i in itertools.count(start=1):
				if not os.path.exists(f'{folder_name}{i:02}/'):
					file_incre = f'{i:02}'
					break
		try:
			Folder_ls = [f'Dock{file_incre}', f'Dock{file_incre}/Receptors', f'Dock{file_incre}/Ligands', 
						f'Dock{file_incre}/output', f'Dock{file_incre}/logs']
			for folder in Folder_ls:
				os.mkdir(folder)
			print('Folders created!!!4')
		except:
			print('Folder exist!!!4')
			pass
	
		return f'Dock{file_incre}'


def prep_input(working_folder:str, pythonsh='', prepLig='', prepReceptor=''):
		
	Receptor_folder= f'{working_folder}/Receptors'
	Lig_folder= f'{working_folder}/Ligands'

	# The prepare_ligand4.py cannot identify paths to ligand. So it is required to change the directory to the coresponding directory.
	os.chdir(Lig_folder)
	# get all ligand PDB files
	Ligs = glob.glob(f'{Lig_folder}/*.pdb')
	for lig in Ligs:
		# Prep Ligand PDBQT files
		os.system(f'{pythonsh} {prepLig} -l {lig}')
		print(f'\n{lig[:-4]}.pdbqt is generated!!!')
	# change the working directory to where the script is
	os.chdir('../../')

	#Prep the receptor PDBQT
	Receps = glob.glob('f{Receptor_folder}/*.pdb')
	for recep in Receps:
		os.system(f'{pythonsh} {prepReceptor} -r {recep} -o {recep[:-4]}.pdbqt')
		print(f'\n{recep[:-4]} is generated!!!!')


def run_vina(working_folder:str, RecpPDBQT:str, LigPDBQT:str, vina='', config='', exhaustiveness=32):

	out_pdbqt = f'{working_folder}/output/{LigPDBQT.rsplit("/",1)[-1]}'
	log = f'{working_folder}/logs/{LigPDBQT.rsplit("/",1)[-1][:-4]}.txt'

	os.system(f'{vina} --receptor {RecpPDBQT} --ligand {LigPDBQT} --config {config} --out {out_pdbqt} --log {log} --exhaustiveness {exhaustiveness}')


#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='''This script is developd to make the work flow of autodock_vina more smooth
				 								Some functions can only be run in PyMOL. ''')
parser.add_argument('-c',
					'--config',
					nargs=1,
					type=str,
					help='Input a configuration file (.txt) that contains the coordinance of the center atom and the box size')
parser.add_argument('-d',
					'--dock', 
					action='store_true',
					help='Initiate docking')
parser.add_argument('-e',
					'--exhaustiveness',
					nargs=1,
					type=int,
					help='Specify an exhaustiveness value, if not provided, 32 is set by default.')
parser.add_argument('-f',
					'--file_incre',
					nargs=1,
					type=int,
					help='Specify a file increment, if not provided, the code will automatically assign one.')
parser.add_argument("-l", 
					'--ligand',
					nargs=1,
					type=str,
					help='Path to one or multiple ligand.pdbqt. If there are multiple ligands, only input the path to Ligand fold')
parser.add_argument('-p',
					'--prep', 
					action='store_true',
					help='Initiate input preparation including preparations for the ligands and receptor')
parser.add_argument("-r", 
					'--receptor',
					nargs=1,
					type=str,
					help='Path to receptor.pdbqt')
parser.add_argument("-s", 
					'--setup_folders',
					nargs='*',
					type=str,
					help='''To set up folders. Ususally it takes two arguments:
							1. folder name
							2. numbering mode
						 ''')
parser.add_argument('-S',
					action='store_true',
					help="For some reason pymol automatically give this flag, without it, script will crash!!!")
parser.add_argument("-w", 
					'--working_directory',
					nargs=1,
					type=str,
					help='Path to the working directory. It is required if -p and -d are not given during folder setup.')

args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------
	
# Here is the main code
def main():

	# Set up paths
	'''
	Check if the running system is WSL or Linux. If it is WSL, by default, the autodock_vina software should be placed on the Desktop.
	If the system is on Linux, autodock_vina should be placed in $HOME (~/). Below section will set the paths of all other softwares
	used in autodock_vina.
	'''
	if os.path.exists('/proc/sys/fs/binfmt_misc/WSLInterop'):
		UserPath = subprocess.check_output('wslpath -u "$(wslvar USERPROFILE)"', shell=True, text=True).strip() # <- take a bit long to run
		UserPath = f'{UserPath}/Desktop'
	else:
		UserPath = f'{os.environ["HOME"]}'

	pythonsh = f'{UserPath}/autodock_vina/mgltools/bin/pythonsh'
	prepLig = f'{UserPath}/autodock_vina/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
	prepReceptor = f'{UserPath}/autodock_vina/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
	vina = f'{UserPath}/autodock_vina/bin/vina'

	# set up working directory
	if args.setup_folders is not None:
		if len(args.setup_folders) == 0:
			if args.file_incre:
				file_incre = f'{args.file_incre[0]:02}'
				working_folder = setup_folders(file_incre=file_incre)
			else:
				working_folder = setup_folders()
		elif len(args.setup_folders) == 1:
			folder_name = args.setup_folders[0]
			if args.file_incre:
				file_incre = f'{args.file_incre[0]:02}'
				working_folder = setup_folders(folder_name, file_incre)
			else:
				working_folder = setup_folders(folder_name)
		elif len(args.setup_folders) == 2:
			folder_name = args.setup_folders[0]
			if args.setup_folders[1] == 'True':
				if args.file_incre:
					file_incre = f'{args.file_incre[0]:02}'
					working_folder = setup_folders(folder_name,file_incre=file_incre, numbering=True)
				else:
					working_folder = setup_folders(folder_name, numbering=True)
			else:
				working_folder = setup_folders(folder_name)

	# Preparing receptor and ligands
	if args.prep:
		if not args.setup_folders:
			if not args.working_directory:
				raise ValueError('Please specify working directory using -w')
			else:
				working_folder = args.working_directory[0]
				prep_input(working_folder, pythonsh=pythonsh, prepLig=prepLig, prepReceptor=prepReceptor)
		else:
			prep_input(working_folder, pythonsh=pythonsh, prepLig=prepLig, prepReceptor=prepReceptor)

	# Run autodock_vina
	if args.dock:
		if not args.receptor:
			raise ValueError('Please input a valid Receptor.pdbqt file')
		elif not args.ligand:
			raise ValueError('Please input a valid Ligand.pdbqt file')
		elif not args.config:
			raise ValueError('''Please input a valid configuration file that contains
			 							the coordinance of the center atom and the box size''')
		else:
			if not args.setup_folders:
				if not args.working_directory:
					raise ValueError('Please specify working directory using -w')
				else:
					working_folder = args.working_directory[0]
					RecepPDBQT = args.receptor[0]
					LigPDBQT = args.ligand[0]
					config_file = args.config[0]
					if args.exhaustiveness:
						exhaustiveness = args.exhaustiveness[0]
						run_vina(working_folder, RecepPDBQT, LigPDBQT, vina, config_file, exhaustiveness)
					else:
						run_vina(working_folder, RecepPDBQT, LigPDBQT, vina, config_file)
			else:
				RecepPDBQT = args.receptor[0]
				LigPDBQT = args.ligand[0]
				config_file = args.config[0]
				if args.exhaustiveness:
					exhaustiveness = args.exhaustiveness[0]
					run_vina(working_folder, RecepPDBQT, LigPDBQT, vina, config_file, exhaustiveness)
				else:
					run_vina(working_folder, RecepPDBQT, LigPDBQT, vina, config_file)

#Below statement prevents the execution of main() when the scrip is run in PyMOL
if __name__ == '__main__' and len(sys.argv) > 1:
	main()
