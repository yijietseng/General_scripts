'''
Created by Yi Jie (Josh) Tseng --- 2/24/2022
Modified by Yi Jie (Josh) Tseng --- 11/11/2023

Description: 
This script is created to prepare receptors, ligands, and to simulate molecular docking using autodock_vina. It can do single, multiple or a library ligand prep and moluecular docking. For more information and instruction, please refer to the README. 
   
'''

import os, itertools, argparse, sys, subprocess, time, threading, re, glob, math, shutil
software_folder = os.getcwd()
home = os.getenv('HOME')

def wall_time(start, end):
    d = end - start
    
    if d < 60:
        d = f'{d} secs'
    elif 60 <= d < 3600:
        d = f'{d / 60} mins'
    elif 3600 <= d < 86400:
        d = f'{d / 3600} hrs'    
    elif 86400 <= d:
        d = f'{d / 86400} days'
    
    return f'[+] Total elapsed time is {d}'
    

def spinner():
    spin = "/-\\|"
    stop_event = threading.Event()

    def _spin():
        while not stop_event.is_set():
            for char in spin:
                sys.stdout.write(f'{char}')
                sys.stdout.flush()
                sys.stdout.write("\b")  # Move the cursor back to overwrite the spin character
                time.sleep(0.1)  # Adjust the delay for the desired speed

    spinner_thread = threading.Thread(target=_spin)
    spinner_thread.start()

    return stop_event


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


def Box(selection, x:float, y:float, z:float, save=None):
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
    #    cmd.delete('Position')
        cmd.delete('Box')
    
    # Get center coordinance. If the input selection is a list, it will automatically use the the coordinaces specified in the list as the center point.
    if isinstance(selection, list):
        pX, pY, pZ = selection
    else:
        pX, pY, pZ = center_coord(selection)
        cmd.delete('Position')
    
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
       4-----------8               ↑
      /|          /|               |
     / |         / |               |
    3-----------7  |               Y	  
    |  |    O   |  |	           |
    |  2--------|--6               o----X----→	
    | /         | /	             ↙
    |/          |/	           Z
    1-----------5            ↙      
    
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


def setup_folders(folder_name=None, file_incre=None, numbering=False):
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
            
            Folder_ls = [f'{folder_name}{file_incre}', f'{folder_name}{file_incre}/Receptors', f'{folder_name}{file_incre}/Ligands', 
                            f'{folder_name}{file_incre}/output', f'{folder_name}{file_incre}/logs']
            for folder in Folder_ls:
                try:
                    os.mkdir(folder)
                    print(f'[+] {folder} created!!!')
                except:
                    print(f'[+] {folder} Folder exist!!!')
                    pass

            return f'{folder_name}{file_incre}'
        
        else:
            if file_incre:
                
                Folder_ls = [f'{folder_name}{file_incre}', f'{folder_name}{file_incre}/Receptors', f'{folder_name}{file_incre}/Ligands', 
                                f'{folder_name}{file_incre}/output', f'{folder_name}{file_incre}/logs']
                
                for folder in Folder_ls:
                    try:
                        os.mkdir(folder)
                        print(f'[+] {folder} created!!!')
                    except:
                        print(f'[+] {folder} Folder exist!!!')  
                        pass
                   
            else:
                Folder_ls = [f'{folder_name}', f'{folder_name}/Receptors', f'{folder_name}/Ligands', 
                             f'{folder_name}/output', f'{folder_name}/logs']
                for folder in Folder_ls:
                    try:
                        os.mkdir(folder)
                        print(f'[+] {folder} created!!!')
                    except:
                        print(f'[+] {folder} Folder exist!!!')  
                        pass

            return folder_name
    else:
        # if no folder_name given, the numbering system is automatically on
        numbering = True
        if file_incre:
            file_incre = file_incre[0]
        else:
            for i in itertools.count(start=1):
                if not os.path.exists(f'Dock{i:02}/'):
                    file_incre = f'{i:02}'
                    break
                
        Folder_ls = [f'Dock{file_incre}', f'Dock{file_incre}/Receptors', f'Dock{file_incre}/Ligands', 
            f'Dock{file_incre}/output', f'Dock{file_incre}/logs']
        for folder in Folder_ls:
            try:
                os.mkdir(folder)
                print(f'[+] {folder} created!!!')
            except:
                print(f'{folder} exists!!!')
                pass
    
        return f'Dock{file_incre}'


def find_max_array_no(working_folder):
    working_folder = working_folder.split('/')[0]
    if glob.glob(f'{working_folder}/Ligands/*.pdbqt'):
        print(len(glob.glob(f'{working_folder}/Ligands/*.pdbqt')))

        return len(glob.glob(f'{working_folder}/Ligands/*.pdbqt'))
    
    else:
        print(len(glob.glob(f'{working_folder}/Ligands/*/*.pdbqt')))

        return len(glob.glob(f'{working_folder}/Ligands/*/*.pdbqt'))
    

def prep_ligand_input_txt(working_folder, all_lig_count):
    print('[+] Prepare ligand input files...')
    all_lig = glob.glob(f'{working_folder}/Ligands/*/*')
    group_size = math.ceil(all_lig_count/5000)
    total_files = math.ceil(all_lig_count/group_size)

    if not os.path.exists(f'{working_folder}/Ligands/Ligand_input_txt'):
        os.mkdir(f'{working_folder}/Ligands/Ligand_input_txt')

    for i in range(1,total_files+1):
        group_lig = all_lig[:group_size]
        with open(f'{working_folder}/Ligands/Ligand_input_txt/LigGrp_{i:04}.txt', 'w') as f:
            f.write('\n'.join(map(str, group_lig)) + '\n')
        all_lig = all_lig[group_size:]
    print('[+] All ligand input files are prepped!')
    

def convert_smi2pdbqt(masterfile):
    #, output='master_lig.pdbqt'
    output = f'{masterfile[:-4]}.pdbqt'
    home_path = os.getenv('HOME')
    print(f'\n[+] Processing {masterfile[:-4]}')
    if os.path.exists(f'{home_path}/fsl_groups/grp_MolecularDock'):
        obabel = f'{home_path}/fsl_groups/grp_MolecularDock/openbabel-2.3.1/build/bin/obabel'
        subprocess.run(f'{obabel} {masterfile} -O {output} --gen3d -h', shell=True, text=True)	
    else:
        subprocess.run(f'obabel {masterfile} -O {output} --gen3d -h', shell=True, text=True)
    print(f'[+] Done processing {masterfile[:-4]}')


def renumber(temp_master):
    count = 0
    with open(f'{temp_master}', 'r') as infile:
        with open(f'{temp_master}2', 'a') as outfile:
            for line in infile:
                if line.startswith('MODEL'):
                    count += 1
                    outfile.write(f'MODEL {count:15}\n')
                else:
                    outfile.write(line)    	
    print('[+] Renumbering complete')

    return count
    

def convert_txt2smi(inputfile):
    '''
    This function converts the .txt file downloaded from the ZINC15 to SMI file.
    '''
    with open(inputfile, 'r') as fin:
        with open(f'{inputfile[:-4]}.smi', 'w') as fout:
            for line in fin:
                name = line.strip().split('\t')[0]
                smi = line.strip().split('\t')[1]
                fout.write(f'{smi} {name}\n')


def split_PDBQT(file=None, working_folder=None, filelimit=None, master='master_lig.pdbqt'):
    
    if not os.path.exists(working_folder):
        os.mkdir(working_folder)
    if not os.path.exists(os.path.join(working_folder,'Ligands/')):
        os.mkdir(os.path.join(working_folder,'Ligands/'))
    if file != None and file.split('.')[-1] == 'pdbqt':
        master_file = file
    else:
        master_file = f'{working_folder}/{master}'

    print('[+] Splitting ligands')
    with open(master_file) as infile:
        count = 0 			# count the number of molecule processed
        in_dir_count = 0 	# count how many in the folder
        dircount = 0 		# count how many folders
        for dircount in itertools.count():
            for line in infile:
                if line.strip().split()[0] == 'MODEL' and line.strip().split()[1] == f'{count+1}':
                    # Make segment folders
                    directory = os.path.join(os.path.join(working_folder,'Ligands/'), f'{dircount+1}')
                    # make subfolders
                    os.makedirs(directory, exist_ok=True)				
                    
                    # create a temp output file and start to crop out molecule info from the master file
                    temp_name = f'temp_{count+1}.pdbqt'
                    out_file = os.path.join(directory, temp_name)
                    with open(out_file, 'w') as outfile:
                        for line in infile:
                            if line.strip() != 'ENDMDL':
                                outfile.write(line)
                                if line.split()[0] == 'REMARK' and line.split()[1] == 'Name':
                                    NewName = os.path.join(directory,f'{line.split()[3]}.pdbqt')
                            else:
                                break
                                
                    os.rename(out_file, NewName)

                    count += 1
                    in_dir_count += 1
                    if in_dir_count >= filelimit:
                        in_dir_count = 0
                        
                        break
            else: 
                break
        
    print('----------\n[+] Done\n')
    return count


def prep_ZINC_lig(filename, working_folder=None, filelimit=None, slurm_gen=False, email=None):
    '''
    To run this function, first go to https://zinc15.docking.org/tranches/home/, click on the 3D tab on the top, and select the tranches of ligands you want to download. Once you have selected the desired tranches, click on the download button, at the bottom of the popup window, select the "AutoDock (*.pdbqt.zt)" for download format and select "WGET" for download method.This function is built to run under Linux environment and to process the ligand in PDBQT format. So when downloading the ZINC15 file, PDBQT and WGET option should be selected. 

    The first argument should be the ZINC15 downloaded file (*.gz.wget) or a processed master file (default as master_lig.pdbqt) with all ligands. The *.gz.wget file contains commands that will allow downloading the ligand pdbqt files. When a file name with a *.gz.wget extension, this function will automatically download the ligands using the commands within the file, it will unzip the downloaded ligands, compile them into a big master file and number each ligand. The master file is for easier file management. 

    After downloading the ligands, if a working folder and a file limit are provided, this function will further segment the master file based on the given file limit. The purposes of segmenting the master file are that 1. Autodock vina only deals with one molecule per PDBQT file 2. lowering the compuational wall time by multi-threading. This function will split the molecules into total#_of_molecules/filelimit groups with filelimit per group. For example, if there are 27 total molecules and a filelimit = 9, it will result in 3 folders and each folder with 9 PDBQT files.
    
    CAUSTION:
    There are A LOT of ligands from https://zinc15.docking.org/tranches/home/. Be mindful of not exceeding the storage and file limit quota. 

    '''
    if not os.path.exists(f'{home}/fsl_groups/grp_MolecularDock'):
        raise SystemError('Please run ZINC Ligand prep on the super computer')

    if '.gz.wget' in filename:
        import gzip

        with open(filename, 'r')as infile:

            for line in infile:
                try:
                    namegz = line.split()[-1]	# get file name with .gz extension, Ex: XX/XXXX/XXXX.YY.pdbqt.gz
                    name = line.split()[-1].split('gz')[0][:-1] # get file name without .gz, Ex: XX/XXXX/XXXX.YY.pdbqt
                    wget = line.strip() # get the wget string in the ZINC15 downloaded file
                    tempfolder = line.split()[-1].split('/',1)[0]

                    cat = f'cat {name} >> temp'
                    
                    # download the pdbqt files
                    os.system(wget)
                
                    # Unzip file
                    with gzip.open(namegz, 'rb') as gzipped_file:
                        with open(name, 'wb') as output:
                            shutil.copyfileobj(gzipped_file, output)
                    print(f'[+] Finish unziping {namegz}')
                    
                    # Add a numbering identifier
                    with open(name) as f:
                        first = f.readline()
                    if first.split()[0] == 'MODEL':
                        os.system(cat)
                    else:
                        os.system('echo "MODEL\t1" >> temp')
                        os.system(cat)
                        os.system('echo "ENDMDL" >> temp')
                except:
                    with open('error', 'a') as e:
                        e.write(line)

                shutil.rmtree(tempfolder)
                print(f'[+] Finish renumbering ligands\n')

        # Fix all the numbers in the temp file
        count = renumber('temp')
        
        # remove all temp files
        os.remove('temp')
        os.rename('temp2', f'{working_folder}/master_lig.pdbqt')
        print('----------\n[+] Ligands are combined into master_lig.pdbqt\n')

    # if a filelimit is given, it will split the molecules into total_molecules/filelimit groups with filelimit per group.
    if filelimit:
        if '.pdbqt' in filename:
            count = split_PDBQT(file=filename, working_folder=working_folder, filelimit=filelimit)
        else:
            split_PDBQT(working_folder=working_folder, filelimit=filelimit)
        
        # Prepare lig file input. This is to reduce the wall time for each slurm job.
        if count >= 5000:
            prep_ligand_input_txt(working_folder, count)

    if slurm_gen:
        print('[+] Writing the slurm submission script')
        print(f'[+] {count} ligands detected')
        update_dock_slurm(email=email, max_array_No=count, working_folder=working_folder)
    
    sys.exit('-----------------------\n[+] Done')


def prep_FDA_lig(input_smi, working_folder=None, filelimit=None, slurm_gen=False, email=None):

    if not os.path.exists(f'{home}/fsl_groups/grp_MolecularDock'):
        raise SystemError('Please run FDA Ligand prep on the super computer')
    
    working_folder = working_folder.split('/')[0]
    #cpuNo = cpu if cpu else os.cpu_count()	
 
    if input_smi.rsplit('.',1)[-1] == 'txt':
        print('[+] Detecting TXT format, converting to SMI format...')
        stop_event = spinner()
        convert_txt2smi(input_smi)
        input_smi = f'{input_smi[:-4]}.smi'
        print('[+] Conversion completed!!!')

  
    if not os.path.exists(working_folder):
        setup_folders(working_folder)
    
 
    os.makedirs('templig', exist_ok=True)

    start = time.time()
    with open(input_smi, 'r') as fin:
        # Get inputfile content
        lines = fin.readlines()
    total_lines = len(lines)
    print(f'[+] {total_lines} ligands detected')
    print('[+] Segmenting ligand files...')

    stop_event = spinner()
    
    file_count = 1
    temp_file_ls = []
    # Creating temp files
    while len(lines) != 0:
        temp_file_name = f'templig/temp{file_count:04}.smi'
        temp_file_ls.append(temp_file_name)
        with open(temp_file_name, 'w') as fout:
            fout.write(lines[0])
        lines = lines[1:]
        file_count += 1
    
    # Set the stop_event to signal the spinner thread to stop
    stop_event.set()
    sys.stdout.write(f'[+] Done segmenting {input_smi} into individual compound smi file\n')
    sys.stdout.flush()
    #print(f'[+] Done segmenting {input_smi} into individual compound smi file')


    # Generating submission script
    update_smi2pdbqt_slurm(max_array_No=len(temp_file_ls))
    
    # submitting updated script using sbatch
    try:
        ID = subprocess.run('sbatch submit_smi2pdbqt.sh', shell=True, capture_output=True, text=True, check=True).stdout.strip()
        slurm_job_ID = int(re.findall(r'\d+', ID)[0])
        print(f'[+] Submitted batch job: {slurm_job_ID}')
    except subprocess.CalledProcessError as e:
        print("Error submitting SLURM job:", e)
    print('[+] Converting SMI to PDBQT...')
    # Sleep for 60 sec before checking completion
    stop_event = spinner()
    
    time.sleep(30)
    task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
    print(f'[+] {len(task)} task(s) remaining...')
    prev_task_length = len(task)
    
    while len(task) > 0:
        time.sleep(20)
        task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
        if len(task) != prev_task_length:
            print(f'[+] {len(task)} task(s) remaining...', flush=True)
            prev_task_length = len(task)  # Update the previous task lengt
        
    stop_event.set()
        
    sys.stdout.write('[+] Done converting all smi files to pdbqt files')
    sys.stdout.flush()
    
    print('\n[+] Distributing PDBQT files...')    
    stop_event = spinner()  
    folder_No = math.ceil(total_lines/filelimit)
    for i in range(1, folder_No+1):
        os.makedirs(f'{working_folder}/Ligands/{i}', exist_ok=True)
            
    folder_count = 1
    in_dir_count = 0
    lig_count = 0
    for temp_pdbqt in glob.glob('templig/*.pdbqt'):
        with open(temp_pdbqt, 'r') as fin:
            for line in fin:
                if line.split()[0] == 'REMARK' and line.split()[1] == 'Name':
                    NewName = os.path.join(working_folder, f'Ligands/{folder_count}/{line.split()[3]}.pdbqt')
                    shutil.move(temp_pdbqt, NewName)
                    lig_count += 1 
                    in_dir_count += 1            
                    if in_dir_count == filelimit:
                        in_dir_count = 0
                        folder_count += 1
                        if folder_count > filelimit:
                            break
                       
    stop_event.set()
    sys.stdout.write(f'[+] Done distributing {lig_count} ligands into {filelimit} folders')
    sys.stdout.flush()
    
    
    print('\n[+] Removing all temporary files...')
    
    shutil.rmtree('templig')

    print(f'[+] Process complete')
    
    if slurm_gen:
        print('\n[+] Updating the slurm submission script')
        if lig_count >= 5000:
            update_dock_slurm(email=email, max_array_No=lig_count, working_folder=working_folder, prep_lig_txt=True)
        else:
            update_dock_slurm(email=email, max_array_No=lig_count, working_folder=working_folder)
    
    sys.exit('-----------------------\n[+] Done')
    

def update_smi2pdbqt_slurm(max_array_No=None, name='submit_smi2pdbqt.sh'):
    with open(name, 'w') as fout:
        fout.write('#!/bin/bash\n\n')
        fout.write('#SBATCH --time=1-00:00:00   # walltime\n')
        fout.write(f'#SBATCH --nodes=1   # number of nodes\n')
        fout.write(f'#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)\n')
        fout.write('#SBATCH --mem-per-cpu=100M   # memory per CPU core\n')
        fout.write(f'#SBATCH -J "smi2pdbqt"   # job name\n')
        fout.write('#SBATCH --output=/dev/null # Do not output slurm ouput\n')
        if max_array_No:
            fout.write(f'#SBATCH --array=1-{max_array_No}\n\n')
        fout.write(f'Lig_array=($(ls templig/*.smi))\n\n')
        fout.write('python3 autodock.py --convert_smi2pdbqt ${Lig_array[$((SLURM_ARRAY_TASK_ID - 1))]}')


def update_dock_slurm(mode='vs', email=None, max_array_No=None, working_folder=None, exhaustiveness=8, name='submit_autodock.sh', prep_lig_txt=False):
    '''
    This function is used to update the slurm submission script. Based on the usage, you can pass different parameters to adjust the script.
    I found that node = 1, cpu = 1, and memory = 100M is sufficient for all applications. 
    
    The submission script will do the following:
    1. Setup all slurm parameters (lines starting with #SBATCH). 
    2. Setup other parameters, such as working folder receptor and configuration file.
    3. Writes submission mechanisms based on the type of applications and the amount of ligands
    
    When submitting the script usually you just need to type: sbatch submit_autodock.sh <path_to_receptor> <path_to_config_file>
    '''
    
    working_folder = f'{working_folder.split("/")[0]}'
    
    with open(name, 'w') as fout:
        
        fout.write('#!/bin/bash --login\n\n'
                   '#SBATCH --time=3-00:00:00   # walltime\n'
                   f'#SBATCH --nodes=1   # number of nodes\n'
                   f'#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)\n')
        
        # Later I found out that 100M usually sufficient, I leave the above code there just in case
        fout.write(f'#SBATCH --mem-per-cpu=1G   # memory per CPU core\n'
                   f'#SBATCH -J "AD_{mode}"   # job name\n')
        if mode == 'vs':
            fout.write(f'#SBATCH --output={working_folder}/SLURM_OUT/slurm-%A_%a.out # Output everything into a temp folder\n')
        if email:
            fout.write(f'#SBATCH --mail-user={email}   # email address\n'
                       '#SBATCH --mail-type=END\n')
        if max_array_No:
            # Calculate group size
            group_size = math.ceil(max_array_No/5000)
            total_group = math.ceil(max_array_No/group_size)
            if max_array_No < 5000:              
                fout.write(f'#SBATCH --array=1-{max_array_No}\n\n')
            elif max_array_No >= 5000:
                fout.write(f'#SBATCH --array=1-{total_group}\n\n')
        
        fout.write(f'Working_folder="{working_folder}"\n'
                   'recep=$1\n'
                   'config=$2\n\n')

        if mode == 'vs':
            if max_array_No < 5000:
                fout.write(f'Lig_array=($(ls {working_folder}/Ligands/*/*))\n\n')
            elif max_array_No >= 5000:
                fout.write('formatted_id=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")\n'
                           'file_path="$Working_folder/Ligands/Ligand_input_txt/LigGrp_$formatted_id.txt"\n\n')
        elif mode == 'exp':
            fout.write(f'Lig_array=($(ls {working_folder}/top_Ligands/*.pdbqt))\n\n')
        else:
            fout.write(f'Lig_array=($(ls {working_folder}/Ligands/*.pdbqt))\n\n')
        
        if mode == 'vs':
            # Check for the existance of SLURM_OUT folder, if not create it.
            fout.write('# Check existance of the temp slurm output folder, if not there, mkdir it.\n'
                    'if [ ! -d "$Working_folder/SLURM_OUT" ]; then\n'
                    '    mkdir -p "$Working_folder/SLURM_OUT"\n'
                    'fi\n\n')
                      
        # If Ligand count < 5000, then use the normal slurm array system
        if max_array_No < 5000:
            fout.write('for i in ${SLURM_ARRAY_TASK_ID};do\n')
            if mode == 'vs':
                fout.write('    ./vina --receptor ${recep} --ligand ${Lig_array[$((SLURM_ARRAY_TASK_ID - 1))]} --out /dev/null --config ${config} --cpu 1 --exhaustiveness 10 | awk -v name="${Lig_array[$((SLURM_ARRAY_TASK_ID - 1))]}" \'$1 == "1" {{print name "\\t" $2; exit}}\'\n ')
            else:
                fout.write('    python3 autodock.py -d -r $recep -l ${Lig_array[$((SLURM_ARRAY_TASK_ID - 1))]} '+'-c $config -e {}\n'.format(exhaustiveness))
            fout.write('done\n\n')
            
        # If Ligand count >= 5000, then use the nested slurm array system                 
        elif max_array_No >= 5000:
            fout.write('while IFS= read -r line; do\n'
                       '''   ./vina --receptor ${recep} --ligand "$line" --out /dev/null --config ${config} --cpu 1 --exhaustiveness 10 | awk -v name="$line" '$1 == "1" {print name "\t" $2; exit}' \n'''
                       'done < "$file_path"\n\n')
            
        # Wait till all subshell processes are finished and concate all temp output into the result file, then remove the temp folder
        #fout.write('# Wait till all subshell processes are finished\n'
        #            'wait\n\n'
        #            '# Concate all results into the actual result file\n'
        #            'cat $Working_folder/SLURM_OUT/slurm-${SLURM_ARRAY_JOB_ID}* >> $Working_folder/logs/Virtural_screen_result.out\n\n')

    if prep_lig_txt:
        prep_ligand_input_txt(working_folder, max_array_No)

    print(f'[+] {name} is updated')
     
     
def prep_input(RecpPDB=None, LigPDB=None):
    
    if not os.path.exists(f'{home}/fsl_groups/grp_MolecularDock'):
        raise SystemError('Please run FDA Ligand prep on the super computer')
    
    # Setting the path to the scripts	
    pythonsh = f'{software_folder}/mgltools/bin/pythonsh'
    prepLig = f'{software_folder}/prepare_ligand4.py'
    prepReceptor = f'{software_folder}/prepare_receptor4.py'
    
  
    if not os.path.exists(pythonsh):
        raise FileNotFoundError(f'{pythonsh} is not found, please modify the path in the script!!!')
    if not os.path.exists(prepLig):
        raise FileNotFoundError(f'{prepLig} is not found, please modify the path in the script!!!')
    if not os.path.exists(prepReceptor):
        raise FileNotFoundError(f'{prepReceptor} is not found, please modify the path in the script!!!')
    
    if RecpPDB:
        working_folder = RecpPDB.rsplit("/",1)[0]
        #Prep the receptor PDBQT if it's not already exists
        if not os.path.exists(f'{RecpPDB.rsplit("/",1)[-1][:-4]}.pdbqt'):
            os.system(f'{pythonsh} {prepReceptor} -r {RecpPDB} -o {working_folder}/{RecpPDB.rsplit("/",1)[-1][:-4]}.pdbqt -A hygdrogens')
            print(f'\n[+]{working_folder}/{RecpPDB.rsplit("/",1)[-1][:-4]}.pdbqt is generated!!!!')

    if LigPDB:
        working_folder = LigPDB.rsplit("/",1)[0]
        Lig_folder = f'{LigPDB.rsplit("/",1)[0]}'
        # The prepare_ligand4.py cannot identify paths to ligand. So it is required to change the directory to the coresponding directory.
        os.chdir(Lig_folder)	

        # Prep ligand pdbqt if it's not already exists
        if not os.path.exists(f'{LigPDB.rsplit("/",1)[-1][:-4]}.pdbqt'):
            os.system(f'{pythonsh} {prepLig} -l {LigPDB} -A hydrogens')
            print(f'\n[+]{LigPDB[:-4]}.pdbqt is generated!!!')
            # change the working directory to where the script is
            os.chdir('../../')


def run_vina(RecepPDBQT:str, LigPDBQT:str=None, config=None, exhaustiveness=32, virtual_screen=False, mode=None):
    import pandas as pd
    
    vina = f'{software_folder}/vina'
    if not os.path.exists(vina):
        raise FileNotFoundError(f'{vina} is not found, please modify the path in the script!!!')

    working_folder = RecepPDBQT.rsplit('/',2)[0]

    if virtual_screen:
        '''
        If virtual_screen is True, it will perform virtual screen on the ZINC15 ligands. To use optimal computational power, cpu is set to 1 and exhaustiveness is set to 8. There will be no option to output structures, instead, it will output an energy file with the top energies for all molecules. From here, you can do a more explicit docking with the best docking molecule by turning off the virtual scrren option and increase the exhaustiveness.
        '''
        start = time.time()
        if not os.path.exists(f'{working_folder}'):
            os.mkdir(f'{working_folder}')
        if not os.path.exists(f'{working_folder}/output'):
            os.mkdir(f'{working_folder}/output')

        # create a result file with header
        with open(f'{working_folder}/logs/Virtural_screen_result.out', 'w') as fout:
            fout.write('Ligand_name\tBinding_Energies_(kcal/mol)\n')

        # submitting updated script using sbatch
        try:
            ID = subprocess.run(f'sbatch submit_autodock.sh {RecepPDBQT} {config}', shell=True, capture_output=True, text=True, check=True).stdout.strip()
            slurm_job_ID = int(re.findall(r'\d+', ID)[0])
            print(f'[+] Submitted batch job: {slurm_job_ID}')
            print('[+] Performing molecular docking...')
        except subprocess.CalledProcessError as e:
            print("Error submitting SLURM job:", e)
            
        stop_event = spinner()
        
        time.sleep(60)
        task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
        print(f'[+] {len(task)} task remaining...')
        
        prev_task_length = len(task)
        
        while len(task) > 0:
            time.sleep(60)
            task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
            if len(task) != prev_task_length:
                print(f'[+] {len(task)} task(s) remaining...', flush=True)
                prev_task_length = len(task)  # Update the previous task lengt
            
        stop_event.set()
        sys.stdout.write('[+] All docking processes are completed')
        sys.stdout.flush()
        
        print('[+] Combining result files...')
        os.system(f'cat {working_folder}/SLURM_OUT/slurm-{slurm_job_ID}_* >> {working_folder}/logs/Virtural_screen_result.out')

        print('[+] Sorting output file...')
        df = pd.read_csv(f'{working_folder}/logs/Virtural_screen_result.out', delimiter='\t')
        sorted_df = df.sort_values(by='Binding_Energies_(kcal/mol)', ascending=True)
        sorted_df.to_csv(f'{working_folder}/logs/Virtural_screen_result.out', sep='\t', index=False)

        print('[+] Removing temporary files...')
        shutil.rmtree(f'{working_folder}/SLURM_OUT')
        
        print('[+] Extracting top scoring ligands...')
        get_tops(working_folder)
        
        stop_event = spinner()
    
        stop_event.set()
        sys.stdout.write('-----------------------------\n[+] Completed')
        sys.stdout.flush()
        try:
            os.remove('core.*')
        except:
            pass
        end = time.time()
    
        wall_time(start, end)

    
    elif mode == 'exp':
        # submitting updated script using sbatch
        print('[+] Initializing docking...')
        try:
            ID = subprocess.run(f'sbatch submit_autodock.sh {RecepPDBQT} {config}', shell=True, capture_output=True, text=True, check=True).stdout.strip()
            slurm_job_ID = int(re.findall(r'\d+', ID)[0])
            print(f'[+] Submitted batch job: {slurm_job_ID}')
            print('[+] Performing explicit molecular docking...')
        except subprocess.CalledProcessError as e:
            print("Error submitting SLURM job:", e)
            
        stop_event = spinner()
        
        time.sleep(15)
        task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
        print(f'[+] {len(task)} task remaining...')
        
        prev_task_length = len(task)
        
        while len(task) > 0:
            time.sleep(15)
            task = subprocess.run(f"squeue -j {slurm_job_ID}", shell=True, capture_output=True, text=True, check=True).stdout.strip().split('\n')[1:]
            if len(task) != prev_task_length:
                print(f'[+] {len(task)} task(s) remaining...', flush=True)
                prev_task_length = len(task)  # Update the previous task length
            
        stop_event.set()
        sys.stdout.write('[+] All docking processes are completed\n')
        sys.stdout.flush()
        
    else:
        # for single docking
        out_pdbqt = f'{working_folder}/output/{LigPDBQT.rsplit("/",1)[-1]}'
        log = f'{working_folder}/logs/{LigPDBQT.rsplit("/",1)[-1][:-6]}.txt'
        
        os.system(f'{vina} --receptor {RecepPDBQT} --ligand {LigPDBQT} --config {config} --out {out_pdbqt} --log {log} --exhaustiveness {exhaustiveness}')


def get_tops(working_folder=None, top=10, result_file=None):
    '''
    If this function is called, it will create a folder and the copy the top ligands into it, default 10 ligands. But the number of top ligands can be specified.
    Additionally, to maintain the system quota at a functional level, it will remove all ligand folders and tar the folder with the top ligands, master_ligand.pdbqt and the result file. And move them into the outFolder. 
    '''
    import pandas as pd
    
    working_folder = working_folder.rsplit('/',1)[0]
    # Specify paths
    if result_file is None:  
        result_path = f'{working_folder}/logs/Virtural_screen_result.out'
    else:
        result_path = result_file
    top_dir = f'{working_folder}/top_Ligands'
    
    if not os.path.exists(top_dir):
        os.mkdir(top_dir)
    
    # copy top ligands
    df = pd.read_csv(result_path, sep='\t')

    for i in range(top):
        lig = df['Ligand_name'][i]
        lig_name = df['Ligand_name'][i].rsplit('/', 1)[-1]
        shutil.copy(lig, f'{top_dir}/{i+1:02}_{lig_name}')
        
    print(f'[+] Top ligands are extracted to the: {top_dir}')
    
    
def package(working_folder=None, tar_out=None, master=None, result_file=None, output_struc=None):
    import getpass
    
    working_folder = working_folder.rsplit('/',1)[0]
    
    username = getpass.getuser()
    outFolder=f'/home/{username}/fsl_groups/grp_MolecularDock/autodock_result' 

    if master is None:
        master_path = f'master_lig.pdbqt'
    else:
        master_path = shutil.move(master, working_folder) 
    if result_file is None:
        result_path = f'logs/'
    else:
        result_path = result_file
    if output_struc is None:
        output_path = f'output/'
    else:
        output_path = output_struc
        
    top_dir = f'top_Ligands/'
        
    # Create output folder
    if not os.path.exists(f'{outFolder}/{working_folder}'):
        os.mkdir(f'{outFolder}/{working_folder}') 
        
    # tar the result package 
    if tar_out is not None:
        original_directory = os.getcwd()
        os.chdir(working_folder)
        
        # Get new paths
        if not os.path.exists(f'{result_path.rsplit("/",1)[0]}'):
            result_path = shutil.move(result_path, working_folder).rsplit('/', 1)[-1]
        if not os.path.exists(f'{output_path.rsplit("/",1)[0]}'):
            output_path = shutil.move(output_path, working_folder).rsplit('/', 1)[-1]
        
        # tar file
        if tar_out.split('.', 1)[-1] == 'tar.gz':
            os.system(f'tar -czf {outFolder}/{working_folder}/{tar_out} {master_path} {top_dir} {result_path} {output_path}')
        else:
            os.system(f'tar -czf {outFolder}/{working_folder}/{tar_out}.tar.gz {master_path} {top_dir} {result_path} {output_path}')
       
        # remove the Ligands and files    
        #shutil.rmtree(f'{working_folder}/Ligands/')
        #shutil.rmtree(f'{top_dir}')
        #os.remove(result_path)
        #os.remove(master_path)
        #os.remove(output_path)
       
       
        os.chdir(original_directory)
        
        try:
            os.mkdir(f'{working_folder}/Ligands')
        except:
            pass
        try:
            os.mkdir(f'{working_folder}/output')
        except:    
            pass
        try:
            os.mkdir(f'{working_folder}/logs')
        except:
            pass
        
        
        
                        
    

#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='''This script is developd to make the work flow of autodock_vina more smooth
                                                 Some functions can only be run in PyMOL. ''')
parser.add_argument('-c',
                    '--config',
                    nargs=1,
                    type=str,
                    metavar='<str>',
                    help='Input a configuration file (.txt) that contains the coordinance of the center atom and the box size')
parser.add_argument('--convert_smi2pdbqt',
                    nargs=1,
                    type=str,
                    metavar='<str>',
                    help='Convert smi file to pdbqt file')
parser.add_argument('-d',
                    '--dock', 
                    action='store_true',
                    help='Initiate docking')
parser.add_argument('-e',
                    '--exhaustiveness',
                    nargs=1,
                    type=int,
                    metavar='<int>',
                    help='Specify an exhaustiveness value, if not provided, 32 is set by default.')
parser.add_argument('--explicit', 
                    action='store_true',
                    help='Turn on explicit docking mode, this is for testing purpose')
parser.add_argument('-f',
                    '--file_incre',
                    nargs=1,
                    type=int,
                    metavar='<int>',
                    help='Specify a file increment, if not provided, the code will automatically assign one.')
parser.add_argument('-F',
                    '--FDA', 
                    nargs="+",
                    help='''Prepare the FDA approved ligands from ZINC15. It uses .smi as input and will output a .pdbqt file with all ligands.
                            It will also segment the ligands into different group for supercomputing. It can also update the submission script 
                            for supercomputing. Argments: input_smi, working_folder=None, filelimit=None, slurm_gen=False, email=None, super=False
                         ''')
parser.add_argument('--get_tops', 
                    nargs="+",
                    help='Get the top ligands. arguments: <working_folder> <# of top structure> <result_file_path>. If result_file is not specified, it will use the default one in the logs folder')
parser.add_argument("-l", 
                    '--ligand',
                    nargs='*',
                    type=str,
                    metavar='<ligand.pdbqt:str>',
                    help='Path to one or multiple ligand.pdbqt. If there are multiple ligands, only input the path to Ligand fold')
parser.add_argument("-m", 
                    '--max_array_no',
                    nargs=1,
                    type=str,
                    metavar='<Working_folder:str>',
                    help='To calculate the maximum number of array necessary for virtual screening')
parser.add_argument('-p',
                    '--prep', 
                    action='store_true',
                    help='Initiate input preparation including preparations for the ligands and receptor')
parser.add_argument('--package', 
                    nargs="+",
                    help='Package the top ligands, master_log.pdbqt and result log file. arguments: <working_folder> <tar_out> <master_file> <result_file> <output_struct>. If the master and the result files are not specified, it will use the default')
parser.add_argument("-r", 
                    '--receptor',
                    nargs=1,
                    type=str,
                    metavar='<receptor.pdbqt:str>',
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
parser.add_argument('-u',
                    '--update_slurm',
                    nargs='+',
                    type=str,
                    help='''Update the slurm submission script.
                            mode:str(vs or norm) <email:str> <max_array_No:int> <working_folder:str> <exhaustiveness=8> <name="submit_autodock.sh> <prep_lig_txt:boolen, default=False>''')
parser.add_argument('-v',
                    '--virtual_screen', 
                    action='store_true',
                    help='Turn on virtual screening mode')
parser.add_argument('-Z',
                    '--ZINC15', 
                    nargs="+",
                    help='''Downloading ZINC15 ligands, unzip, combine, renumber and split the argument should be the file name
                             of the file downloaded from ZINC15 (*.gz.wget) or the combined master file with .pdbqt extension. The
                            second and third arguments are a working folder and a file limit value. Both are optional. They are 
                            only used when splitting the master file is needed.
                            <filename:str> <working_folder:str> <filelimit:int> <slurm_gen:boolen> <email:str>''')

args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------
    
# Here is the main code
def main():

    # set up working directory
    if args.setup_folders is not None:
        if len(args.setup_folders) == 0:
            if args.file_incre:
                file_incre = f'{args.file_incre[0]:02}'
                setup_folders(file_incre=file_incre)
            else:
                setup_folders()
        elif len(args.setup_folders) == 1:
            folder_name = args.setup_folders[0]
            if args.file_incre:
                file_incre = f'{args.file_incre[0]:02}'
                setup_folders(folder_name, file_incre)
            else:
                setup_folders(folder_name)
        elif len(args.setup_folders) == 2:
            folder_name = args.setup_folders[0]
            if args.setup_folders[1] == 'True':
                if args.file_incre:
                    file_incre = f'{args.file_incre[0]:02}'
                    setup_folders(folder_name,file_incre=file_incre, numbering=True)
                else:
                    setup_folders(folder_name, numbering=True)
            else:
                setup_folders(folder_name)

    # Calculate maximum array no
    if args.max_array_no:
        find_max_array_no(args.max_array_no[0])
    
    # Prep ZINC15 ligands
    if args.ZINC15:
        if len(args.ZINC15) == 1:
            prep_ZINC_lig(filename=args.ZINC15[0])
        elif 1 < len(args.ZINC15) < 4:
            prep_ZINC_lig(filename=args.ZINC15[0],working_folder=args.ZINC15[1],filelimit=int(args.ZINC15[2]))
        elif len(args.ZINC15) == 4:
            prep_ZINC_lig(filename=args.ZINC15[0],working_folder=args.ZINC15[1],filelimit=int(args.ZINC15[2]),slurm_gen=eval(args.ZINC15[3]))
        elif len(args.ZINC15) > 4:
            prep_ZINC_lig(filename=args.ZINC15[0],working_folder=args.ZINC15[1],filelimit=int(args.ZINC15[2]),slurm_gen=eval(args.ZINC15[3]), email=args.ZINC15[4])
            
    if args.FDA:
        if len(args.FDA) == 1:
            prep_FDA_lig(input_smi=args.FDA[0])
        elif 1 < len(args.FDA) == 3:
            prep_FDA_lig(input_smi=args.FDA[0], working_folder=args.FDA[1], filelimit=int(args.FDA[2]))
        elif len(args.FDA) == 5:
            prep_FDA_lig(input_smi=args.FDA[0], working_folder=args.FDA[1], filelimit=int(args.FDA[2]), slurm_gen=eval(args.FDA[3]), email=args.FDA[4])

    if args.convert_smi2pdbqt:
        convert_smi2pdbqt(args.convert_smi2pdbqt[0]) 
            
    if args.update_slurm:
        if len(args.update_slurm) == 6:
            update_dock_slurm(mode=args.update_slurm[0], email=args.update_slurm[1], max_array_No=int(args.update_slurm[2]), working_folder=args.update_slurm[3], 
                              exhaustiveness=int(args.update_slurm[4]), name=args.update_slurm[5])
        elif len(args.update_slurm) == 5:
            update_dock_slurm(mode=args.update_slurm[0], email=args.update_slurm[1], max_array_No=int(args.update_slurm[2]), working_folder=args.update_slurm[3], 
                              exhaustiveness=int(args.update_slurm[4]))
        elif len(args.update_slurm) == 4:
            update_dock_slurm(mode=args.update_slurm[0], email=args.update_slurm[1], max_array_No=int(args.update_slurm[2]), working_folder=args.update_slurm[3])


    # Preparing receptor and ligands
    if args.prep:
        if args.receptor:
            if args.ligand:
                if len(args.ligand) > 1:
                    for lig in args.ligand:
                        prep_input(RecpPDB=args.receptor[0], LigPDB=lig)
                else:
                    prep_input(RecpPDB=args.receptor[0], LigPDB=args.ligand[0])
            else:
                prep_input(RecpPDB=args.receptor[0])
        else:
            if args.ligand:
                if len(args.ligand) > 1:
                    for lig in args.ligand:
                        prep_input(LigPDB=lig)
                else:
                    prep_input(LigPDB=args.ligand[0])
            else:
                raise ValueError('Please input a valid receptor or ligand PDB files!!!')
    
    # Run autodock_vina
    if args.dock:
        if not args.receptor:
            raise ValueError('Please input a valid Receptor.pdbqt file')
        elif not args.ligand and not args.virtual_screen:
            raise ValueError('Please input a valid Ligand.pdbqt file')
        elif not args.config:
            raise ValueError('''Please input a valid configuration file that contains
                                         the coordinance of the center atom and the box size''')
        else:
            if args.virtual_screen:
                # get email
                with open('submit_autodock.sh', 'r') as fin:
                    for line in fin:
                        if 'email' in line:
                            email = line.split('=')[-1].split('   # ')[0]
                            
                run_vina(RecepPDBQT=args.receptor[0], config=args.config[0], virtual_screen=True)
                
                # After the run, do another explicit docking with the top ligands
                working_folder = args.receptor.rsplit('/',2)[0]
                
                top_ligs = glob.glob(f'{working_folder}/top_Ligands/*')
                
                update_dock_slurm(mode='exp', email=email ,max_array_No=len(top_ligs), working_folder=working_folder, exhaustiveness=256)
                if not os.path.exists(f'{working_folder}/output'):
                    os.mkdir(f'{working_folder}/output')
                run_vina(RecepPDBQT=args.receptor[0], config=args.config[0], exhaustiveness=args.exhaustiveness[0], mode='exp')
            else:
                if args.exhaustiveness:
                    run_vina(RecepPDBQT=args.receptor[0], LigPDBQT=args.ligand[0], config=args.config[0], exhaustiveness=args.exhaustiveness[0])
                else:
                    run_vina(RecepPDBQT=args.receptor[0], LigPDBQT=args.ligand[0], config=args.config[0])
    
    if args.explicit:
        
        email = input('Please input your email: ').strip()
        # After the run, do another explicit docking with the top ligands
        working_folder = args.receptor[0].rsplit('/',2)[0]
        if not os.path.exists(f'{working_folder}/top_Ligands'):
            raise FileNotFoundError(f'{working_folder}/top_Ligands not found')
        
        top_ligs = glob.glob(f'{working_folder}/top_Ligands/*')
        update_dock_slurm(mode='exp', email=email ,max_array_No=len(top_ligs), working_folder=working_folder, exhaustiveness=10)
        run_vina(RecepPDBQT=args.receptor[0], config=args.config[0], exhaustiveness=args.exhaustiveness[0], mode='exp')
    
    if args.get_tops:
        if len(args.get_tops) == 1:
            get_tops(working_folder=args.get_tops[0])
        elif len(args.get_tops) == 2:
            get_tops(working_folder=args.get_tops[0], top=int(args.get_tops[1]))
        else:
            get_tops(working_folder=args.get_tops[0], top=int(args.get_tops[1]), result_file=args.get_tops[2])
            
    if args.package:
        if len(args.package) < 2:
            raise ValueError('Please enter <working_folder>, <tar_out>')
        elif len(args.package) == 2:
            package(args.package[0], args.package[1])
        elif len(args.package) == 3:
            package(args.package[0], args.package[1], args.package[2])
        elif len(args.package) == 4:
            package(args.package[0], args.package[1], args.package[2], args.package[3])
        else:
            package(args.package[0], args.package[1], args.package[2], args.package[3], args.package[4])
            
            
#Below statement prevents the execution of main() when the scrip is run in PyMOL
if __name__ == '__main__' and len(sys.argv) > 1:
    main()


