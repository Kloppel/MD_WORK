#class NAMD also entails CHARMM executables
class NAMD():
    def create_executables(model_text, filename, ending):
        #open file as writeable
        with open(f'{filename}.{ending}', 'w') as f:
            #write
            f.write(model_text)
        return

    def execute_charmm(molecule_name, path_text):
        #create individual bash script name
        exec_name = f'build_{molecule_name}.inp'
        #create bash script text with shebang and charmm execution command
        bash_text = f'#!/bin/bash \n \n cd {path_text} charmm -i {exec_name} -o {molecule_name}.log \n'
        #write bash file executing charmm commands to a file
        NAMD.create_executables(model_text=bash_text, filename=f'bash_charmm_{molecule_name}', ending="sh")
        #execute charmm file through bash script
        subprocess.call(f'bash_charmm_{molecule_name}.sh')
        return

    def execute_namd(molecule_name, **kwargs):
        # kwargs: heating_type, run, nodes, 
        #create individual bash script name
        exec_name = f'build_{molecule_name}.inp'

        #for heating runs
        if kwargs.get('heating_type') is not None:
            path_text = organizer.create_folders_heat(heating_type=heating_type)
            bash_text = f'#!/bin/bash \n \n cd {path_text} \n/home/pbuser/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +idlepoll +p{nodes} +devices 0 {exec_name} > {molecule_name}.log \n cd ../ \n'
            #write bash file executing charmm commands to a file
            NAMD.create_executables(model_text=bash_text, filename=f'bash_namd_{molecule_name}', ending="sh")
            #execute charmm file through bash script, bash script will move into specified results folder depending on if it is a heating or production run
            subprocess.call(f'bash_namd_{molecule_name}.sh')

        #for production runs
        if kwargs.get('run') is not None:
            if kwargs.get('nodes') is not None:
                for run in runs:
                    path_text = organizer.create_folders_prod(run=run)
                    bash_text = f'#!/bin/bash \n \n cd {path_text} \n/home/pbuser/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +idlepoll +p{nodes} +devices 0 {exec_name} > {molecule_name}.log \n cd ../ \n'
                    #write bash file executing charmm commands to a file
                    NAMD.create_executables(model_text=bash_text, filename=f'bash_namd_{molecule_name}', ending="sh")
                    #execute charmm file through bash script, bash script will move into specified results folder depending on if it is a heating or production run, then execute the calculation and move out of it
                    subprocess.call(f'bash_namd_{molecule_name}.sh')
            else:
                sys.exit('You said we were executing a production run, but did not specify how many nodes to use. \nSince running a production run on a CUDA on only one node (which otherwise I would have to assume) \nis kinda senseless, the input of "nodes" is necessary.')
        if kwargs.get() is None:
            raise TypeError('Please give either a heating_type (vel_resc or constr_relax) or a number of runs for a production run. \nThis function needs either one of those to work, otherwise the output will be None.')
        return
    
    def build_molecule(molecule_name, redox_state, **kwargs):
        if molecule_name != "cco":
            print('Molecule is not cco, no prior data is found. Please consider that you have to specify variables \n protein_parameter_file, protein_coordinates, membrane_coordinates, water_coordinates \n since none of these coordinates are deposited by the author.')
        #define standard topology files if none are given
        if kwargs.get('clean_topology_file') is None :
            clean_topology_file = "\"top_alw_clean.inp\""
        if kwargs.get('lipid_topology_file') is None:
            lipid_topology_file = "\"top_all36_lipid.rtf\""
        #define standard parameter files if none are given
        if kwargs.get('protein_parameter_file') is None and molecule_name=="cco":
            protein_parameter_file = "\"par_all22_prot_plus_heme_and_Cu.inp\""
        if kwargs.get('lipid_parameter_file') is None:
            lipid_parameter_file = "\"par_all36_lipid.prm\""
        #define cco standard coordinates
        if kwargs.get('protein_coordinates') is None and molecule_name=="cco":
            protein_coordinates = f'read sequence coor resid name "KetteA" \nGENERATE  ACHAIN SETUP \nread coor card name "KetteA" \nread sequence coor resid name "KetteB" \nGENERATE BCHAIN SETUP \nread coor card name "KetteB" \nread sequence coor resid name "HemA3" \nGENERATE EHEMEA3 SETUP \nread coor card name "HemA3" \nread sequence coor resid name "HemA" \nGENERATE GHEMEA SETUP \nread coor card name "HemA" \nread sequence coor resid name "Metalle" \nGENERATE METAL SETUP \nread coor card name "Metalle" \nread sequence coor resid name "Hydroxyl" \nGENERATE FEOH SETUP \nread coor card name "Hydroxyl" \nread sequence coor resid name "waterc" \nGENERATE HOHCU SETUP \nread coor card name "waterc" \n'
        #define standard membrane coordinates for cco
        if kwargs.get('membrane_coordinates') is None and molecule_name=="cco":
            membrane_coordinates = f'read sequence coor resid name "membrane" \nGENERATE MEMBRANE SETUP \nread coor card name "membrane" \n'
        #define standard waters coordinates for cco
        if kwargs.get('water_coordinates') is None and molecule_name=="cco":
            water_coordinates = f'read sequence coor resid name "H2OA"\nGENERATE PAH2O SETUP noangle nodihedral\nread coor card name "H2OA"\nread sequence coor resid name "H2OB"\nGENERATE QBH2O SETUP noangle nodihedral\nread coor card name "H2OB"\nread sequence coor resid name "waterm"\nGENERATE SWAT SETUP noangle nodihedral\nread coor card name "waterm"\nread sequence coor resid name "waterbox"\nGENERATE WAT SETUP noangle nodihedral\nread coor card name "waterbox"\n'
        #CcO redox states in order of reaction in a dictionary, so that dict values are inserted as text-form 
        #patches into build file.
        cco_redox_state_collection= {
            "Pm": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBP2 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H4 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISO  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n' ,
            "Pr": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBPN METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H4 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISO  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "P->F": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBP2 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A343 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISO  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "F": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBF2 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H4 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISO  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "F->O": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBF2 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A343 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISO  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "O": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBF2 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H3 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISE  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "O->E": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBT4 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H3 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISE  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUB2 METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n' , 
            "E": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CB1T METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H3 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISE  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUB2 METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "E->R": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CBT4 METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3W2 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISW  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n', 
            "R": f'!Patches\n!Sulfur bridge\nPATCH DISU ACHAIN 64 ACHAIN 88 SETUP\n!Protonation patches\nPATCH GLUP ACHAIN 286 SETUP\nPATCH ASPP ACHAIN 407 SETUP\n! PATCH LSN ACHAIN 362 SETUP ! depending if Lys from K-channel should be deprotonated\n! PATCH GLUP BCHAIN 101 SETUP ! depending if Glu101 at K-channel entrance should be protonated\n!generate angles and dihedrals\nAUTOGENERATE ANGLES DIHEDRALS\n!charge patches set up the charge for the redox-active cofactors -> here patches used for O->E transition\nPATCH AHE3 GHEMEA 2  ACHAIN 102 ACHAIN  421 SETUP ! heme a\nPATCH CA21 METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP ! copper A\nPATCH CB1T METAL 1 ACHAIN 333 ACHAIN 334 HOHCU 1 ACHAIN 288 ACHAIN 284 SETUP ! copper B\nPATCH A3H3 EHEMEA3 2 ACHAIN 419 FEOH 1 SETUP ! heme a3\n!bonds towards heme and copper \nPATCH PHEM  ACHAIN 419  EHEMEA3 2 SETUP\nPATCH PHEM  ACHAIN  102 GHEMEA 2 SETUP\nPATCH PHE2  ACHAIN 421 GHEMEA 2 SETUP\nPATCH EISW  FEOH 1 EHEMEA3 2 SETUP\nPATCH CUBP METAL 1 ACHAIN 284 ACHAIN 333 ACHAIN 334 HOHCU 1 SETUP\nPATCH CUAP METAL 2 METAL 3 BCHAIN 217 BCHAIN 252 BCHAIN 254 BCHAIN 256 BCHAIN 260 BCHAIN 263 SETUP\n'
        }
        if molecule_name == "cco":
            redox_state_patch = cco_redox_state_collection[redox_state]
        if molecule_name == "cco" and redox_state not in ("Pm", "Pr", "P->F", "F", "F->O", "O", "O->E", "E", "E->R", "R"):
            raise ValueError("\n \nWARNING \nUnknown Redox State of CcO given, viable states are: \nPm, Pr, P->F, F, F->O, O, O->E, E, E->R, R\n WARNING \n \n")
        model_text = f'* THIS IS THE CHARMM PROGRAM THAT BUILDS {molecule_name}\n* \n \n!read topology files params\nOPEN READ UNIT 42 CARD NAME {clean_topology_file} !"top_alw_clean.inp"\nREAD rtf CARD UNIT 42\nclose unit 42\n \nOPEN READ UNIT 42 CARD NAME {lipid_topology_file} !"top_all36_lipid.rtf"\nREAD rtf CARD UNIT 42 APPEND\nclose unit 42\n \nOPEN READ UNIT 42 CARD NAME {protein_parameter_file} !"par_all22_prot_plus_heme_and_Cu.inp"\nREAD PARAMETERS UNIT 42 CARD\nclose unit 42\n \nOPEN READ UNIT 42 CARD NAME {lipid_parameter_file} !"par_all36_lipid.prm"\nREAD PARAMETERS UNIT 42 CARD APPEND\nclose unit 42\n \n!read coordinates\n{protein_coordinates}\n{membrane_coordinates}\n \n!make hydrogen\nhbuild\n \n{redox_state_patch}\n \n!read water coordinates\n{water_coordinates}\n \n! Another hbuild to find possible H-atoms that may be included since the patches \nhbuild\n!!!!!!!MINIMIZATION OF HYDROGENS AFTER THE STATE CHANGE WOULD BE RECOMMENDED!!!!!!!\n \nOPEN WRITE UNIT 08 CARD NAME "{molecule_name}.crd" !"cco_and_water.crd"\nWRITE COORDINATES UNIT 08 CARD\nCLOSE UNIT 08\n \nOPEN WRITE UNIT 08 CARD NAME "{molecule_name}.crd" !"cco_and_water.psf"\nWRITE PSF UNIT 08 CARD\nCLOS UNIT 08\n \nOPEN WRITE UNIT 08 CARD NAME "{molecule_name}.crd" !"cco_and_water.psf.xplor"\nWRITE PSF XPLOR UNIT 08 CARD\nCLOS UNIT 08\n \nOPEN WRITE UNIT 08 CARD NAME "{molecule_name}.crd" !"cco_and_water.pdb"\nWRITE COORDINATES UNIT 08 PDB \nCLOS UNIT 08 \n \nstop\n'
        #Now create the executable from this created text.
        path_text = NAMD.create_executables(model_text=model_text, filename=f'build_{molecule_name}', ending="inp")
        #NAMD.execute_charmm(molecule_name, path_text)
        return 
    

    def velocity_rescaling_heating(molecule_name, num_steps, timestep, output_name, output_energies, temperature, **kwargs):
        # kwargs: _1_4scaling, switch_dist, cutoff, pair_list_dist, margin, rigid_tolerance, rigid_iterations, langevin, langevin_damping, langevin_temp, fixed_atoms_file, pme, pme_x, pme_y, pme_z, scaling_high, scaling_low
        if int(num_steps) % 2 == 0:
            pass
        else:
            raise ValueError("Error: EvenIntegerTotalStepsError: Total number of steps must be an even integer, since the number of steps must be of integer size, \nand even because constraint scaling will happen twice with different constraints. \nThe number of steps is calculated from the input for num_steps (total steps), divided by two.\nConverting to the closest even integer now.")
            num_steps = math.ceil(num_steps/2.)*2
        #Standard CcO Parameter and Topology Files
        if kwargs.get('lipid_eq') is None and molecule_name == "cco":
            lipid_eq = "2gsm_md_membrane"
        if kwargs.get("protein_parameter_file") is None and molecule_name == "cco":
            protein_parameter_file = "par_all22_prot_plus_heme_and_Cu"
        if kwargs.get("lipid_parameter_file") is None and molecule_name == "cco":
            lipid_parameter_file = "par_all36_lipid"
        if kwargs.get('constraints_file') is None and molecule_name == "cco":
            constraints_file = "namd_cons_for_heat"
        if kwargs.get("fixed_atoms_file") is None and molecule_name == "cco":
            fixed_atoms_file = "../namd_heat.pdb"
        if kwargs.get('constraints_file') is None and molecule_name == "cco":
            constraints_file = "namd_cons_for_heat"
        
        #Standard Values from old CcO Simulations
        #Force Field Params
        if kwargs.get('dielectric') is None and molecule_name == "cco":
            dielectric = "1.0"
        if kwargs.get('_1_4scaling') is None:
            _1_4scaling = "1.0"
        if kwargs.get('cutoff') is None:
            cutoff = "12.0"
        if kwargs.get('switch_dist') is None:
            switch_dist = "10.0"
        if kwargs.get('pair_list_dist') is None:
            pair_list_dist = "13.5"
        if kwargs.get('margin') is None:
            margin = "1.0"
        if kwargs.get('rigid_tolerance') is None:
            rigid_tolerance = '0.00001'
        if kwargs.get('rigid_iterations') is None:
            rigid_iterations = "100"
        
        #Thermostat and Barostat Standard Values
        if kwargs.get("langevin") is None:
            langevin = "on"
        if kwargs.get("langevin_damping") is None:
            langevin_damping = "1"
        if kwargs.get("pme") is None:
            pme = "yes"
        if kwargs.get('margin_on') is None:
            margin_on = "#"

        #Output Standard Values
        if kwargs.get("output_pressure") is None:
            output_pressure = output_energies
        if kwargs.get('DCD_freq') is None:
            DCD_freq = output_energies
        if kwargs.get("XST_freq") is None:
            XST_freq = output_energies

        #Constraint scaling 
        if kwargs.get("scaling_high") is None:
            scaling_high = "0.70"
        if kwargs.get("scaling_low") is None:
            scaling_low = "0.50"

        #Velocity Rescaling standard Parameters
        if kwargs.get("t0") is None:
            t0 = "0"
        if kwargs.get('min_steps') is None:
            min_steps = "500"
        steps = str(0.5*int(num_steps))

        #Create Restraint pdb from vmd scripting language .tcl
        if kwargs.get('selection1') is None:
            selection1 = 'name CA C O N'
        if kwargs.get('selection2') is None:
            selection2 = 'all and not (hydrogen or ion or water or (name CA C O N))'
        if kwargs.get('pdb_name') is None:
            pdb_name = constraints_file
        vmd_text = f'# load into vmd \n#vmd scripting data type : .tcl \n# your data \nmol load psf cco_and_water.psf pdb cco_and_water.pdb \n#HEAT: \nset all [atomselect top all] \n#segname hemes metals festhalten, alle anderen constrained \n#set fix [atomselect top "(resname AUR S5CD S5CP S6OH NIC F3S F4S SF4)"] \n#$fix set occupancy 1 \nset res1 [atomselect top "{selection1}"] \nset res2 [atomselect top "{selection2}"] \n$all set beta 0 \n$all set occupancy 0 \n$res1 set beta 1.00 \n$res2 set beta 0.50 \n$all writepdb {pdb_name}.pdb \n#quit \n#Befehöl zum Ausführen \n#vmd -dispdev text -e restrain.tcl >restrain.out'
        NAMD.create_executables(model_text=vmd_text, filename='restrain', ending='tcl')
        os.system('vmd -dispdev text -e restrain.tcl >restrain.out')
        #Actual Model
        model_text = f'#############################################################\n## JOB DESCRIPTION                                         ##\n#############################################################\n\n# Minimization and Equilibration of\n# CcO in a membrane in a Water Box\n\n\n#############################################################\n## ADJUSTABLE PARAMETERS                                   ##\n#############################################################\n\nstructure         {molecule_name}.psf.xplor\ncoordinates       {molecule_name}.pdb\nset outputname     {molecule_name}_heated\n\nfirsttimestep      {t0}\nextendedSystem   {lipid_eq}.restart.xsc # for newstart from lipid-equilibration\n\n#############################################################\n## SIMULATION PARAMETERS                                   ##\n#############################################################\n\n# Input\nparaTypeCharmm	    on\nparameters          {protein_parameter_file}.inp \nparameters          {lipid_parameter_file}.prm \n\n# Force-Field Parameters\ndielectric	    {dielectric};# Value of the dielectric constant (added)\nexclude             scaled1-4   ;# Exclude/scale local (along the sequence)\n1-4scaling          {_1_4scaling};# Scale factor for (i,i+4) EL interactions\ncutoff              {cutoff};# In PME, cutoff dictates the separation between long and short range forces and doesnt simply turn off interactions. \n switching           on		;# Energy Switch VdW interactions and partition EL into local and non-local contributions\nswitchdist          {switch_dist};# distance at which to activate switching function for van der Waals (and electrostatic in sph. cut.)\npairlistdist        {pair_list_dist};# CUTNB in charmm\n\n\n# Integrator Parameters\ntimestep            {timestep};# \n#nonbondedFreq       1	;# timesteps between nonbonded evaluation. Positive integer factor of fullElectFrequency.\nfullElectFrequency  2	;# number of timesteps between full electrostatic evaluations Acceptable Values. positive integer factor of stepspercycle\n#margin             0.0 ;# Extra distance used in selecting patches\nstepspercycle       10	;# INBFRQ in charmm, Frequency of updating Verlet list (in integration steps) Def. 20\n### SETTLE\nrigidBonds all 		;# Use SHAKE on solute (and SETTLE on water according to useSettle)\nuseSettle on		;#Apply SETTLE (implemented for water)\nrigidTolerance {rigid_tolerance}	;# Desired accuracy in maintaining SHAKEed bond lengths\nrigidIterations {rigid_iterations};# Maximum number of SHAKE iterations\n\n\n\n\n### Particle Mesh Ewald\nPME on ;# Use PME for electrostatic calculation - USe only factors of 2,3,5\n#PMEGridSizeX 100\n#PMEGridSizeY 100\n#PMEGridSizeZ 100\nPMETolerance 0.000001;\nPMEInterpOrder 4 ;\nPMEGridSpacing      1.0\n\n\n\n\n# Constant Temperature Control\nlangevin            on    ;# do langevin dynamics\nlangevinDamping     1     ;# damping coefficient (gamma) of 5/ps\nlangevinTemp        {temperature}\nlangevinHydrogen    off    ;# dont couple langevin bath to hydrogens\nwrapWater on  	# Are water molecules translated back to the unit cell (purely cosmetic option, has no effect on simulations\nwrapAll on     # Are other molecules such as proteins translated back\n{margin_on}margin 5 # enlarge error acceptance for membrane\n\n# Constant Pressure Control (variable volume)\nuseGroupPressure      yes ;# needed for rigidBonds\nuseFlexibleCell       yes  # changed for membrane\nuseConstantRatio       yes  # changed for membrane\n\nlangevinPiston        {langevin} #on\nlangevinPistonTarget  1.01325 ;#  in bar -> 1 atm\nlangevinPistonPeriod  200.  # changed from 100 for membrane\nlangevinPistonDecay   50.\nlangevinPistonTemp    {temperature}\n\n\n# Output\noutputName          $outputname\n\nrestartfreq         {num_steps};\ndcdfreq             {DCD_freq}\nxstFreq             {XST_freq}\noutputEnergies      {output_energies}\noutputPressure      {output_pressure}\n\n\n#############################################################\n## EXTRA PARAMETERS                                        ##\n#############################################################\n\n#############################################################\n## HARMONIC CONSTRAINTS                                     ##\n#############################################################\n\n#halte alle Atome constrained außer Hs\nconstraints     on\nconsRef         {constraints_file}.pdb #namd_cons_for_heat.pdb\nconsKFile       {constraints_file}.pdb #namd_cons_for_heat.pdb\nconsKCol        B\n\n#############################################################\n## EXECUTION SCRIPT                                        ##\n#############################################################\n\n# Minimization\nminimize            {min_steps} #500\nreinitvels          {temperature}\n\n\n# Heating (total 50 ps)\n#set IHTFRQ 500\n#set TEMINC 6\n\nVELOCITY RESCALING\nfor {{set i 6}} {{$i <= $temperature}} {{incr i $TEMINC }} {{\n    langevinTemp $i\n    reinitvels $i\n    run $IHTFRQ\n}}\n\n# START\n# Langevin dynamics (heating)\n#constraintScaling \n'
        NAMD.create_executables(model_text=model_text, filename=f'rescale_velocities_{molecule_name}', ending="inp")
        #NAMD.execute_namd()
        return

    def constraint_relaxation_heating(molecule_name, num_steps, timestep, output_name, output_energies, temperature, run_high, run_low, **kwargs):
        # kwargs: _1_4scaling, switch_dist, cutoff, pair_list_dist, margin, rigid_tolerance, rigid_iterations, langevin, langevin_damping, langevin_temp, fixed_atoms_file, pme, pme_x, pme_y, pme_z, scaling_high, scaling_low
        if int(num_steps) % 2 == 0:
            pass
        else:
            raise ValueError("Error: EvenIntegerTotalStepsError: Total number of steps must be an even integer, since the number of steps must be of integer size, \nand even because constraint scaling will happen twice with different constraints. \nThe number of steps is calculated from the input for num_steps (total steps), divided by two.\nConverting to the closest even integer now.")
            num_steps = math.ceil(num_steps/2.)*2
        #Standard CcO Parameter and Topology Files
        if kwargs.get('lipid_eq') is None and molecule_name == "cco":
            lipid_eq = "2gsm_md_membrane"
        if kwargs.get("protein_parameter_file") is None and molecule_name == "cco":
            protein_parameter_file = "par_all22_prot_plus_heme_and_Cu"
        if kwargs.get("lipid_parameter_file") is None and molecule_name == "cco":
            lipid_parameter_file = "par_all36_lipid"
        if kwargs.get('constraints_file') is None and molecule_name == "cco":
            constraints_file = "namd_cons_for_heat"
        if kwargs.get("fixed_atoms_file") is None and molecule_name == "cco":
            fixed_atoms_file = "../namd_heat.pdb"
        if kwargs.get('constraints_file') is None and molecule_name == "cco":
            constraints_file = "namd_cons_for_heat"
        
        #Standard Values from old CcO Simulations
        #Force Field Params
        if kwargs.get('dielectric') is None and molecule_name == "cco":
            dielectric = "1.0"
        if kwargs.get('_1_4scaling') is None:
            _1_4scaling = "1.0"
        if kwargs.get('cutoff') is None:
            cutoff = "12.0"
        if kwargs.get('switch_dist') is None:
            switch_dist = "10.0"
        if kwargs.get('pair_list_dist') is None:
            pair_list_dist = "13.5"
        if kwargs.get('margin') is None:
            margin = "1.0"
        if kwargs.get('rigid_tolerance') is None:
            rigid_tolerance = '0.00001'
        if kwargs.get('rigid_iterations') is None:
            rigid_iterations = "100"
        
        #Thermostat and Barostat Standard Values
        if kwargs.get("langevin") is None:
            langevin = "on"
        if kwargs.get("langevin_damping") is None:
            langevin_damping = "1"
        if kwargs.get("pme") is None:
            pme = "yes"
        if kwargs.get('margin_on') is None:
            margin_on = "#"

        #Output Standard Values
        if kwargs.get("output_pressure") is None:
            output_pressure = output_energies
        if kwargs.get('DCD_freq') is None:
            DCD_freq = output_energies
        if kwargs.get("XST_freq") is None:
            XST_freq = output_energies

        #Constraint scaling 
        if kwargs.get("scaling_high") is None:
            scaling_high = "0.70"
        if kwargs.get("scaling_low") is None:
            scaling_low = "0.50"

        #Constraint Scaling standard Parameters
        if kwargs.get("t0") is None:
            t0 = "0"
        if kwargs.get('min_steps') is None:
            min_steps = "500"
        
        steps = str(0.5*int(num_steps))
        model_text = f'#############################################################\n## JOB DESCRIPTION                                         ##\n#############################################################\n\n# Minimization and Equilibration of\n# CcO in a membrane in a Water Box\n\n\n#############################################################\n## ADJUSTABLE PARAMETERS                                   ##\n#############################################################\n\nstructure         {molecule_name}.psf.xplor\ncoordinates       {molecule_name}.pdb\nset outputname     {molecule_name}_heated\n\nfirsttimestep      {t0}\nextendedSystem   {lipid_eq}.restart.xsc # for newstart from lipid-equilibration\n\n#############################################################\n## SIMULATION PARAMETERS                                   ##\n#############################################################\n\n# Input\nparaTypeCharmm	    on\nparameters          {protein_parameter_file}.inp \nparameters          {lipid_parameter_file}.prm \n\n# Force-Field Parameters\ndielectric	    {dielectric};# Value of the dielectric constant (added)\nexclude             scaled1-4   ;# Exclude/scale local (along the sequence)\n1-4scaling          {_1_4scaling};# Scale factor for (i,i+4) EL interactions\ncutoff              {cutoff};# In PME, cutoff dictates the separation between long and short range forces and doesnt simply turn off interactions. \n switching           on		;# Energy Switch VdW interactions and partition EL into local and non-local contributions\nswitchdist          {switch_dist};# distance at which to activate switching function for van der Waals (and electrostatic in sph. cut.)\npairlistdist        {pair_list_dist};# CUTNB in charmm\n\n\n# Integrator Parameters\ntimestep            {timestep};# \n#nonbondedFreq       1	;# timesteps between nonbonded evaluation. Positive integer factor of fullElectFrequency.\nfullElectFrequency  2	;# number of timesteps between full electrostatic evaluations Acceptable Values. positive integer factor of stepspercycle\n#margin             0.0 ;# Extra distance used in selecting patches\nstepspercycle       10	;# INBFRQ in charmm, Frequency of updating Verlet list (in integration steps) Def. 20\n### SETTLE\nrigidBonds all 		;# Use SHAKE on solute (and SETTLE on water according to useSettle)\nuseSettle on		;#Apply SETTLE (implemented for water)\nrigidTolerance {rigid_tolerance}	;# Desired accuracy in maintaining SHAKEed bond lengths\nrigidIterations {rigid_iterations};# Maximum number of SHAKE iterations\n\n\n\n\n### Particle Mesh Ewald\nPME on ;# Use PME for electrostatic calculation - USe only factors of 2,3,5\n#PMEGridSizeX 100\n#PMEGridSizeY 100\n#PMEGridSizeZ 100\nPMETolerance 0.000001;\nPMEInterpOrder 4 ;\nPMEGridSpacing      1.0\n\n\n\n\n# Constant Temperature Control\nlangevin            on    ;# do langevin dynamics\nlangevinDamping     1     ;# damping coefficient (gamma) of 5/ps\nlangevinTemp        {temperature}\nlangevinHydrogen    off    ;# dont couple langevin bath to hydrogens\nwrapWater on  	# Are water molecules translated back to the unit cell (purely cosmetic option, has no effect on simulations\nwrapAll on     # Are other molecules such as proteins translated back\n{margin_on}margin 5 # enlarge error acceptance for membrane\n\n# Constant Pressure Control (variable volume)\nuseGroupPressure      yes ;# needed for rigidBonds\nuseFlexibleCell       yes  # changed for membrane\nuseConstantRatio       yes  # changed for membrane\n\nlangevinPiston        {langevin} #on\nlangevinPistonTarget  1.01325 ;#  in bar -> 1 atm\nlangevinPistonPeriod  200.  # changed from 100 for membrane\nlangevinPistonDecay   50.\nlangevinPistonTemp    {temperature}\n\n\n# Output\noutputName          $outputname\n\nrestartfreq         {num_steps};\ndcdfreq             {DCD_freq}\nxstFreq             {XST_freq}\noutputEnergies      {output_energies}\noutputPressure      {output_pressure}\n\n\n#############################################################\n## EXTRA PARAMETERS                                        ##\n#############################################################\n\n#############################################################\n## HARMONIC CONSTRAINTS                                     ##\n#############################################################\n\n#halte alle Atome constrained außer Hs\nconstraints     on\nconsRef         {constraints_file}.pdb #namd_cons_for_heat.pdb\nconsKFile       {constraints_file}.pdb #namd_cons_for_heat.pdb\nconsKCol        B\n\n#############################################################\n## EXECUTION SCRIPT                                        ##\n#############################################################\n\n# Minimization\nminimize            {min_steps} #500\nreinitvels          {temperature}\n\n\n# Heating (total 50 ps)\n#set IHTFRQ 500\n#set TEMINC 6\n\n\n# START\n# Langevin dynamics (heating)\nconstraintScaling {scaling_high}\nrun {steps}\n\n# Langevin dynamics (heating)\nconstraintScaling {scaling_low}\nrun {steps}\n'
        NAMD.create_executables(model_text=model_text, filename=f'rescale_velocities_{molecule_name}', ending="inp")
        #NAMD.execute_namd()
        return

    def production_run(molecule_name, steps_per_run, runs, timestep, output_name, output_energies, temperature, **kwargs):
        tot_steps = int(steps_per_run) * int(runs)
        if tot_steps % 2 == 0:
            pass
        else:
            raise ValueError("Error: EvenIntegerTotalStepsError: Total number of steps must be an even integer, since the number of steps must be of integer size, \nand even because constraint scaling will happen twice with different constraints. \nThe number of steps is calculated from the input for num_steps (total steps), divided by two.\nConverting to the closest even integer now.")
            num_steps = math.ceil(num_steps/2.)*2
        #maybe read number of usable nodes from the system somehow?
        if kwargs.get('nodes') is None:
            nodes="20"

        #Standard CcO Parameter and Topology Files
        if kwargs.get('lipid_eq') is None and molecule_name == "cco":
            lipid_eq = "2gsm_md_membrane"
        if kwargs.get("protein_parameter_file") is None and molecule_name == "cco":
            protein_parameter_file = "par_all22_prot_plus_heme_and_Cu"
        if kwargs.get("lipid_parameter_file") is None and molecule_name == "cco":
            lipid_parameter_file = "par_all36_lipid"
        if kwargs.get('constraints_file') is None and molecule_name == "cco":
            constraints_file = "namd_cons_for_heat"
        if kwargs.get("fixed_atoms_file") is None and molecule_name == "cco":
            fixed_atoms_file = "../namd_heat.pdb"

        #Standard Values from old CcO Simulations
        #Force Field Params
        if kwargs.get('dielectric') is None and molecule_name == "cco":
            dielectric = "1.0"
        if kwargs.get('_1_4scaling') is None:
            _1_4scaling = "1.0"
        if kwargs.get('cutoff') is None:
            cutoff = "12.0"
        if kwargs.get('switch_dist') is None:
            switch_dist = "10.0"
        if kwargs.get('pair_list_dist') is None:
            pair_list_dist = "13.5"
        if kwargs.get('margin') is None:
            margin = "1.0"
        if kwargs.get('rigid_tolerance') is None:
            rigid_tolerance = '0.00001'
        if kwargs.get('rigid_iterations') is None:
            rigid_iterations = "100"

        #Thermostat and Barostat Standard Values
        if kwargs.get("langevin") is None:
            langevin = "on"
        if kwargs.get("langevin_damping") is None:
            langevin_damping = "1"
        if kwargs.get("pme") is None:
            pme = "yes"
        if kwargs.get('margin_on') is None:
            margin_on = "#"
        
        #Output Standard Values
        if kwargs.get("output_pressure") is None:
            output_pressure = output_energies
        if kwargs.get('DCD_freq') is None:
            DCD_freq = output_energies
        if kwargs.get("XST_freq") is None:
            XST_freq = output_energies

        #EXECUTION with a loop that generates the execution script, moves into the execution folder, 
        #then executes the script, moves out of the execution folder, and adds 1 to run, to indicate
        #a finished run. 
        run = 0
        while run <= int(runs):
            t0 = (run * int(steps_per_run))
            organizer.create_folders_prod(run=run)
            model_text = f'#############################################################\n#############################################################\n######################### DESCRIPTION #######################\n#############################################################\n#############################################################\n# Production Run of\n# {molecule_name} in a Water Box\n#############################################################\n######             ADJUSTABLE PARAMETERS               ######\n#############################################################\nstructure         ../{molecule_name}.psf.xplor\ncoordinates       ../{molecule_name}.pdb\nbinvelocities   ../{molecule_name}_heated.vel\nbincoordinates  ../{molecule_name}_heated.coor\n# remove temperature $temperature command below\n# remove manual perodic boundary definition\n# and margin!\n# and timesteps! -> 5 more ns\nset temperature    {temperature}\nset outputname     {molecule_name}_prod\nfirsttimestep      {t0}\nextendedSystem  ../{molecule_name}_heated.xsc  # for newstart from lipid-equilibration\n#############################################################\n######               SIMULATION PARAMETERS           ########\n#############################################################\n# Input\nparaTypeCharmm	    on\nparameters          ../{protein_parameter_file}.inp\nparameters          ../{lipid_parameter_file}.prm\n#############################################################\n######               FORCE FIELD PARAMS              ########\n#############################################################\ndielectric	    {dielectric}		;# Value of the dielectric constant (added)\nexclude             scaled1-4   ;# Exclude/scale local (along the sequence)\n1-4scaling          {_1_4scaling}		;# Scale factor for (i,i+4) EL interactions\ncutoff              {cutoff}		;# In PME, cutoff dictates the separation between long and short range forces and doesnt simply turn off interactions.\nswitching           on		;# Energy Switch VdW interactions and partition EL into local and non-local contributions\nswitchdist          {switch_dist}		;# distance at which to activate switching function for van der Waals (and electrostatic in sph. cut.)\npairlistdist        {pair_list_dist}	;# CUTNB in charmm\n#############################################################\n######                INTEGRATOR PARAMS                ######\n#############################################################\ntimestep            {timestep}	;\n#nonbondedFreq       1	;# timesteps between nonbonded evaluation. Positive integer factor of fullElectFrequency.\nfullElectFrequency  2	;# number of timesteps between full electrostatic evaluations Acceptable Values. \n#positive integer factor of stepspercycle\n#margin             {margin} ;# Extra distance used in selecting patches\nstepspercycle       10	;\n# INBFRQ in charmm, Frequency of updating Verlet list (in integration steps) Def. 20\n#############################################################\n######                  SETTLE/SHAKE                   ######\n#############################################################\nrigidBonds all 		;# Use SHAKE on solute (and SETTLE on water according to useSettle)\nuseSettle on		;#Apply SETTLE (implemented for water)\nrigidTolerance {rigid_tolerance}	;# Desired accuracy in maintaining SHAKEed bond lengths\nrigidIterations {rigid_iterations}	;# Maximum number of SHAKE iterations\n#############################################################\n######               PARTICLE MESH EWALD               ######\n#############################################################\n#PME off\nPME on ;\n# Use PME for electrostatic calculation - USe only factors of 2,3,5\n#PMEGridSizeX 100\n#PMEGridSizeY 100\n#PMEGridSizeZ 100\n##ischrgd=1, ;\n# NAMD doesnt force neutralization of char_charge\n#PMETolerance 0.000001 ;# Def 1E-6\n#PMEInterpOrder 4 ;# Def 4 (cubic+1)\nPMEGridSpacing      1.0\n#############################################################\n######              CONSTANT TEMP CONTROL              ######\n#############################################################\nlangevin            {langevin}    ;# do langevin dynamics\nlangevinDamping     {langevin_damping}     ;# damping coefficient (gamma) of 5/ps\nlangevinTemp        {temperature}\nlangevinHydrogen    off    ;# dont couple langevin bath to hydrogens\n#############################################################\n######          PERIODIC BOUNDARY CONDITIONS           ######\n#############################################################\n# this block defines periodic boundary conditions\n#cellBasisVector1  100.0    0.     0.	\n# Direction of the x basis vector for a unit cell\n#cellBasisVector2    0.   100.0    0.	\n# Direction of the y basis vector for a unit cell\n#cellBasisVector3    0.     0.   118.4	\n# Direction of the z basis vector for a unit cell\n#cellOrigin          0.     0.     0.	\n# Position of the unit cell center\nwrapWater on  \n# Are water molecules translated back to the unit cell (purely cosmetic option, has no effect on simulations\nwrapAll on    # Are other molecules such as proteins translated back\n{margin_on}margin 5 # enlarge error acceptance for membrane\n#############################################################\n######            CONSTANT PRESSURE CONTROL            ######\n######                (VARIABLE VOLUME)                ######\n#############################################################\nuseGroupPressure      yes ;# needed for rigidBonds\nuseFlexibleCell       yes  # changed for membrane\nuseConstantRatio       yes  # changed for membrane\nlangevinPiston        on\nlangevinPistonTarget  1.01325 ;\n#  in bar -> 1 atm\nlangevinPistonPeriod  200.  \n# changed from 100 for membrane\nlangevinPistonDecay   50.\nlangevinPistonTemp    $temperature\n#############################################################\n######                                                 ######\n######                     OUTPUT                      ######\n######                                                 ######\n#############################################################\noutputName          $outputname\nrestartfreq         50000 ;\n#500steps = every 1ps\ndcdfreq             {DCD_freq}\nxstFreq             {XST_freq}\noutputEnergies      {output_energies}\noutputPressure      {output_pressure}\n#############################################################\n######               EXTRA PARAMETERS                  ######\n#############################################################\n\n#############################################################\n######                                                 ######\n######                    RUN                          ######\n######                                                 ######\n#############################################################\n#Description: Specifies whether or not harmonic constraints are\n#active. If it is set to off, then no harmonic constraints are\n#computed.\n#If it is set to on, then harmonic constraints are calculated using\n#the values specified by the parameters consref, conskfile, conskcol,\n#and consexp.\n#constraints on\n#Description: PDB file to use for reference positions for harmonic\n#constraints. Each atom that has an active constraint will be\n#constrained\n#about the position specified in this file.\n#consref /scratch/scratch/awoelke/md_cco/membrane/solvate/cco_3hb3_in_water.pdb\n#Description: PDB file to use for force constants for harmonic constraints.\n#conskfile /scratch/scratch/awoelke/md_cco/membrane/solvate/min_constrains.pdb\n#Description: Column of the PDB file to use for the harmonic\n#constraint force constant. This parameter may specify any of the\n#floating point\n#fields of the PDB file, either X, Y, Z, occupancy, or beta-coupling\n#(temperature-coupling). Regardless of which column is used, a value of\n#0\n#indicates that the atom should not be constrained. Otherwise, the\n#value specified is used as the force constant for that atoms\n#restraining potential\n#Acceptable Values: X, Y, Z, O, or B\n#conskcol O\n#Description: The harmonic constraint energy function is multiplied\n#by this parameter, making it possible to gradually turn off\n#constraints during\n#equilibration. This parameter is used only if constraints is set to on.\n#constraintScaling 20.0\n#############################################################\n######                   FIX ATOMS                     ######\n#############################################################\n#Description: Specifies whether or not fixed atoms are present.\n#fixedAtoms on\n#Description: Specifies whether or not forces between fixed atoms are calculated. This option is required to turn fixed atoms off in the middle\n#of a simulation. These forces will affect the pressure calculation, and you should leave this option off when using constant pressure if the\n#coordinates of the fixed atoms have not been minimized. The use of constant pressure with significant numbers of fixed atoms is not recommended.\n#fixedAtomsForces off\n#Description: PDB file to use for the fixed atom flags for each atom. If this parameter is not specified, then the PDB file specified by coordinates\n#is used.\n#fixedAtomsFile /scratch/scratch/awoelke/md_cco/membrane/solvate/min_constrains.pdb\n#Description: Column of the PDB file to use for the containing fixed atom parameters for each atom. The coefficients can be read from any floating point\n#column of the PDB file. A value of 0 indicates that the atom is not fixed.\n#Acceptable Values: X, Y, Z, O, or B\n#Default Value: O\n#fixedAtomsCol O\n#############################################################\n######               EXECUTION SCRIPT                  ######\n#############################################################\n#START\n run {steps_per_run} ; ' 
            NAMD.create_executables(model_text=model_text, filename=f'prod{run}/{molecule_name}_prod_run{run}', ending="inp")
            NAMD.execute_namd(molecule_name="cco", run=run, nodes=nodes)
            run += 1
        return

