class organizer():
    def create_folders_prod(run):
        path_text = f'prod{run}/'
        if not os.path.exists(path_text):
            os.makedirs(path_text)
        return path_text
    
    def create_folders_heat(heating_type):
        path_text = f'heat_{heating_type}/'
        if not os.path.exists(path_text):
            os.makedirs(path_text)
        return path_text

    def create_result_folders(**kwargs):
        if not os.path.exists(f'res/logs/'):
            os.makedirs(f'res/logs/')
        if not os.path.exists(f'res/restart_files/'):
            os.makedirs(f'res/restart_files/')
        return

    def move_files(run, mol_name):
        #move regular files into /res, the results folder
        endings = ["coor", "inp", "dcd", "vel", "xsc", "xst"]
        for ending in endings:
            try:
                dummy = Path(f'prod{run}/{mol_name}.{ending}')
                dummy.resolve(strict=True)
            except FileNotFoundError:
                print(f'did not find prod{run}/{mol_name}.{ending}. Breaking Execution now, \nsince a necessary file has not been created by NAMD.')
                break
            else:
                print(f'prod{run}/{mol_name}.{ending} was found, moving files and moving on.')
                shutil.move(f'prod{run}/{mol_name}.{ending}', f'res/{mol_name}_{run}.{ending}')

        #move logs into /log instead of regular folder
        log_endings = ["log"]
        for ending in log_endings:
            try:
                dummy = Path(f'prod{run}/{mol_name}.{ending}')
                dummy.resolve(strict=True)
            except FileNotFoundError:
                print(f'did not find prod{run}/{mol_name}.{ending}. Breaking Execution now, \nsince a necessary file has not been created by NAMD.')
                break
            else:
                print(f'prod{run}/{mol_name}.{ending} was found, moving files and moving on.')
                shutil.move(f'prod{run}/{mol_name}.{ending}', f'res/{mol_name}_{run}.{ending}')

        #move restart files into /restart_files, the restarting files folder
        restart_endings= ["restart.coor", "restart.vel", "restart.xsc"]
        for ending in restart_endings:
            try:
                dummy = Path(f'prod{run}/{mol_name}.{ending}')
                dummy.resolve(strict=True)
            except FileNotFoundError:
                print(f'did not find prod{run}/{mol_name}.{ending}. Breaking Execution now, \nsince a necessary file has not been created by NAMD.')
                break
            else:
                print(f'prod{run}/{mol_name}.{ending} was found, moving files and moving on.')
                shutil.move(f'prod{run}/{mol_name}.{ending}', f'res/{mol_name}_{run}.{ending}')
        return
    
    #def concatenate(steps_per_run, number_of_runs):
    #    #call bash script concatenating all files
    #    subprocess.call("concat.sh")
    #    return
#
    #def delete_folders(run):
    #    if os.path.exists(f'prod{run}/'):
    #        os.removedirs(f'prod{run}/')
    #    return
    