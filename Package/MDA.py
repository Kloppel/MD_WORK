import kb2plus_package
class MDA():
    def empty_universe(n_atoms):
        #check if n_atoms is a correct datatype
        if n_atoms is None:
            raise TypeError
        if isinstance(n_atoms, int) is False:
            raise TypeError
        #create empty universe with a maximum number of atoms n_atoms
        universe = mda.Universe.empty(n_atoms=n_atoms, trajectory=True)
        #desired output, when calling universe object:
        #<Universe with (n_atoms) atoms>
        return universe
    
    def create_waterverse(waters, **kwargs):
        print('filling universe with water.')
        n_residues = waters
        n_atoms = n_residues * 3
        resindices = np.repeat(range(n_residues), 3)
        assert len(resindices) == n_atoms
        segindices = [0] * n_residues
        waterverse = MDA.empty_universe(n_atoms=n_atoms)
        waterverse.add_TopologyAttr('name', ['O', 'H1', 'H2']*n_residues)
        waterverse.add_TopologyAttr('type', ['O', 'H', 'H']*n_residues)
        waterverse.add_TopologyAttr('resname', ['SOL']*n_residues)
        waterverse.add_TopologyAttr('resid', list(range(1, n_residues+1)))
        waterverse.add_TopologyAttr('segid', ['SOL'])
        #positions can be simply assigned:
        h2o = np.array([[ 0,        0,       0      ],  # oxygen
                        [ 0.95908, -0.02691, 0.03231],  # hydrogen
                        [-0.28004, -0.58767, 0.70556]]) # hydrogen
        grid_size, spacing, coordinates = 10, 8, []
        # translating h2o coordinates around a grid
        for i in range(n_residues):
            x = spacing * (i % grid_size)
            y = spacing * ((i // grid_size) % grid_size)
            z = spacing * (i // (grid_size * grid_size))
            xyz = np.array([x, y, z])
            coordinates.extend(h2o + xyz.T)
        coord_array = np.array(coordinates)
        assert coord_array.shape == (n_atoms, 3)
        waterverse.atoms.positions = coord_array
        assert not hasattr(sol, 'bonds')
        bonds = []
        for o in range(0, n_atoms, 3):
            bonds.extend([(o, o+1), (o, o+2)])
        waterverse.add_TopologyAttr('bonds', bonds)
        return waterverse

    def add_segment(universe, segid, atom_selection):
        seg = universe.add_Segment(segID = segid)
        seg_atoms = universe.select_atoms(atom_selection)
        seg.residues.segments = segment
        return universe

    def tile_universe(universe, n_x, n_y, n_z):
        box = universe.dimensions[:3]
        copied = []
        for x in range(n_x):
            for y in range(n_y):
                for z in range(n_z):
                    u_ = universe.copy()
                    move_by = box*(x, y, z)
                    u_.atoms.translate(move_by)
                    copied.append(u_.atoms)

        new_universe = mda.Merge(*copied)
        new_box = box*(n_x, n_y, n_z)
        new_universe.dimensions = list(new_box) + [90]*3
        return new_universe

    def import_structure(structure, waterverse, **kwargs):
        #import protein
        protiverse = mda.Universe(structure)
        #center around origin
        water_cog_ = waterverse.atoms.center_of_geometry()
        waterverse.atoms.positions -= water_cog_
        water_cog = waterverse.atoms.center_of_geometry()
        protein_cog_ = protiverse.atoms.center_of_geometry()
        protiverse.atoms.positions -= protein_cog_
        protein_cog = protiverse.atoms.center_of_geometry()
        print('Water Center of geometry '+ water_cog)
        print('Protein Center of geometry '+ protein_cog)
        for pos in water_cog:
            if pos > 0.1:
                print('Aligning origins resulted in unusually large values. Please consider looking at the center of geometry operations of the waterverse.')
        for pos in protein_cog:
            if pos > 0.1:
                print('Aligning origins resulted in unusually large values. Please consider looking at the center of geometry operations of the protiverse.')
        #create merged universe
        universe = mda.Merge(waterverse.atoms, protiverse.atoms)
        #find the size of our universe as 1.5 times the largest size in one direction
        dist = 1.5*(max(u.atoms.positions[:, 0]) - min(u.atoms.positions[:, 0]))
        #fix the size of our universe, 90 for 90° angles (box-shaped universe)
        universe.dimensions = [dist, dist, dist, 90, 90, 90]
        #tile out the universe (explicit PBCs)
        if kwargs.get('tile') is not None:
            if kwargs.get('n_x') is not None:
                if kwargs.get('n_y') is not None:
                    if kwargs.get('n_z') is not None:
                        universe = MDA.tile_universe(universe=universe, n_x=n_x, n_y=n_y, n_z=n_z)
            else:
                sys.exit('Tiling out the system was requested, but no specifications for tiling have been given. Please enter n_x, n_y, n_z values as function arguments.')
        return universe

    def align_structures(universe):
        #on-the-fly version of aligning protein structures, so that memory is saved.
        #takes a universe object, the object in MDAnalysis in which all simulated/created molecules live, the space, 
        #and time if it's a time-series file
        transforms = [trans.unwrap(),
                    trans.center_in_box(),
                    trans.wrap()]
        universe.trajectory.add_transformations(*transforms)
        return universe

    def rmsd_alignment(universe1, universe2, selection, match_by_mass):
        merged_universe = mda.Merge(universe1.atoms, universe2.atoms)
        #match_by_mass will match atoms by mass if True, and match atoms by number if False
        rmsds = align.alignto(universe1, universe2, select=selection, match_atoms=match_by_mass)
        if match_by_mass == True:
            print(f'Automatically aligning universe1 and universe2 has found the RMSD of ' + rmsds + f'\nThe alignment was carried out through matching atom masses.')
        if match_by_mass ==False:
            print(f'Automatically aligning universe1 and universe2 has found the RMSD of ' + rmsds + f'\nThe alignment was carried out through matching atom numbering.')
        return rmsds, merged_universe
    
    def rmsf_alignment(universe1, universe2, selection, output_file, match_by_mass):
        merged_universe = mda.Merge(universe1.atoms, universe2.atoms)
        align.AlignTraj(universe1, universe2, select=selection, filename=output_file, match_atoms=match_by_mass).run()
        return merged_universe

    def rmsd(universe1, universe2, structure_trajectory, molecule_name, **kwargs):
        #kwargs: atom_selection_struc1, atom_selection_struc2, atom_selection_traj, ref_frame, weights
        if kwargs.get('ref_frame') is None: 
            ref_frame = 0
        if kwargs.get('atom_selection_struc1') is None:
            atom_selection_struc1 = 'all'
        if kwargs.get('atom_selection_struc2') is None:
            atom_selection_struc2 = 'all'
        if kwargs.get('atom_selection_traj') is None:
            atom_selection_traj = 'all'
        if kwargs.get('weights') is None:
            weights = 'mass'

        #for comparing structures
        if structure_trajectory == False:
            RMS = rms.rmsd(universe1.select_atoms(atom_selection_struc1).positions,
                universe2.select_atoms(atom_selection_struc2).positions,
                center=True,
                superposition=True)
            print('Comparing two structures gives the RMSD ' + RMS)
            df = pd.DataFrame(RMS.rmsd, columns=['Frame', 'Time [ns]', atom_selection_struc1, atom_selection_struc2])
            ax = df.plot(x='Frame', y=[atom_selection_struc1, atom_selection_struc2, 'Backbone'])
            ax.set_ylabel(r'RMSD [$\AA$]')
            ax.savefig(f'RMSD_{molecule_name}_{atom_selection_struc1}_{atom_selection_struc2}.png')
        #for comparing trajectory and structure
        if structure_trajectory == True:
            RMS = rms.RMSD(universe1, universe2, select=atom_selection_traj, groupselections=[], ref_frame=ref_frame, weights=weights)
            RMS.run()
            print('Comparing trajectory and structure gives mass-wigthed RMSD ' + RMS)
            df = pd.DataFrame(RMS.rmsd, columns=['Frame', 'Time [ns]', atom_selection_traj])
            ax = df.plot(x='Frame', y=[atom_selection_traj, 'Backbone'])
            ax.set_ylabel(r'RMSD [$\AA$]')
            ax.savefig(f'RMSD_{molecule_name}_{atom_selection_traj}_weigthed_mass.png')
        return
    
    def pairwise_rmsd_self(universe, molecule_name, **kwargs):
        if kwargs.get('selection') is None:
            selection = 'name CA'
        aligner = MDA.align_structures(universe=universe)
        matrix = diffusionmap.DistanceMatrix(universe, select={selection}).run()
        plt.imshow(matrix.dist_matrix, cmap='viridis')
        plt.xlabel('Frame')
        plt.ylabel('Frame')
        plt.colorbar(label=r'RMSD ($\AA$)')
        plt.savefig(f'{molecule_name}_pairwise_self_RMSD_{selection}.png')
        return matrix

    def pairwise_rmsd(universe1, universe2, molecule_name, **kwargs):
        if kwargs.get('selection') is None:
            selection = 'name CA'
        prmsd = np.zeros((len(universe1.trajectory), len(universe2.trajectory)))
        for i, frame_open in enumerate(universe1.trajectory):
            r = rms.RMSD(universe2, universe1, select=selection,
                 ref_frame=i).run()
            prmsd[i] = r.rmsd[:, -1]  # select 3rd column with RMSD values
        plt.imshow(prmsd, cmap='viridis')
        plt.xlabel('Frame (universe2)')
        plt.ylabel('Frame (universe1)')
        plt.colorbar(label=r'RMSD ($\AA$)')
        plt.savefig(f'{molecule_name}_pairwise_RMSD_{selection}.png')
        return prmsd

    def rmsf_trajectory(universe1, universe2, name_u1, name_u2, molecule_name):
        #u = universe1
        #ref = universe2
        average = align.AverageStructure(universe1, universe1, select='protein and name CA', ref_frame=0).run()
        ref = average.universe
        aligner = align.AlignTraj(universe1, universe2, select='protein and name CA', in_memory=False).run()
        c_alphas = universe1.select_atoms('protein and name CA')
        RMS = rms.RMSF(c_alphas).run()
        plt.plot(c_alphas.resids, RMS.rmsf)
        plt.xlabel('Residue number')
        plt.ylabel(r'RMSF ($\AA$)')
        plt.legend()
        plt.savefig(f'{molecule_name}_RMSF_traj_CA.png')
        return RMS

    def distance_between_groups(universe1, universe2, molecule_name, **kwargs):
        #box dimensions: [10, 10, 10, 90, 90, 90]
        #kwargs: selection, PBC, PBC->box, array, selection_u1, selection_u2, atom-wise
        print('Function calculating and plotting the distance between groups from two universes, or one universe if specified for both input arguments \nuniverse1 and universe2. If a specific selection of non-equal atom lists is desired, the keyword array \n has to be specified.')
        if kwargs.get('atom-wise') is None:
            if kwargs.get('selection') is None:
                selection = 'name CA'
                print('Assuming to calculate distances between C-alphas of the two universes, since no selection has been specified.')
            if kwargs.get('PBC') is False:
                group1 = universe1.select_atoms(selection)
                group2 = universe2.select_atoms(selection)
                resids1, resids2, dist = distances.dist(ca1, ca2, offset=0)
            if kwargs.get('PBC') is True:
                if kwargs.get('box') is None:
                    box=universe1.dimensions
                    resids1, resids2, dist_box = distances.dist(ca1, ca2, box=box)
                else: 
                    resids1, resids2, dist_box = distances.dist(ca1, ca2, box=box)
            plt.plot(resids1, dist)
            plt.ylabel(r'CA distance ($\AA$)')
            plt.legend()
            plt.savefig(f'{molecule_name}_group_distance_{selection}.png')
            if kwargs.get('array') is True:
                if kwargs.get('selection_u1') is not None:
                    if kwargs.get('selection_u2') is not None:
                        selection1 = universe1.select_atoms(selection_u1)
                        selection2 = universe2.select_atoms(selection_u2)
                        if kwargs.get('center_of_mass') is True:
                            selection1 = selection1.center_of_mass(compound='residues')
                            selection2 = selection2.center_of_mass(compound='residues')
                        atoms_in_selection1 = len(selection1)
                        atoms_in_selection2 = len(selection2)
                    else: 
                        raise TypeError('SelectionError: selection_u2 is necessary.')
                else:
                    raise TypeError('SelectionError: selection_u1 is necessary.')
            dist_arr = distances.distance_array(selection1.positions, selection2.positions, box=u.dimensions)
            fig, ax = plt.subplots()
            im = ax.imshow(dist_arr, origin='upper')
            tick_interval = 5
            ax.set_yticks(np.arange(atoms_in_selection1)[::tick_interval])
            ax.set_xticks(np.arange(atoms_in_selection2)[::tick_interval])
            ax.set_yticklabels(selection1.resids[::tick_interval])
            ax.set_xticklabels(selection2.resids[::tick_interval])
            plt.ylabel(selection_u1)
            plt.xlabel(selection_u2)
            plt.title(f'Distance between groups: {selection_u1} and {selection_u2}')
            cbar = fig.colorbar(im)
            cbar.ax.set_ylabel('Distance (Angstrom)')
            plt.savefig(f'{molecule_name}_group_distance_array.png')
        if kwargs.get('atom-wise') is not None:
            if kwargs.get('center_of_mass') is None:
                atoms = universe1.select_atoms(selection1)
            if kwargs.get('center_of_mass') is not None:
                atoms = universe1.atoms.center_of_mass(compound='residues')
            number_of_atoms = len(atoms)
            self_distances = distances.self_distance_array(atoms.positions)
            self_distances_array = np.zeros((number_of_atoms, number_of_atoms))
            index = np.triu_indices_from(self_distances_array)
            self_distances_array[index] = self_distances
            self_distances_array.T[index] = self_distances
            fig, ax = plt.subplots()
            im = ax.pcolor(atoms.resids, atoms.resids, self_distances_array)
            ax.set_aspect('equal')
            plt.ylabel('Residue IDs')
            plt.xlabel('Residue IDs')
            plt.title(f'Distance between all atoms in {molecule_name}')
            cbar = fig.colorbar(im)
            cbar.ax.set_ylabel(r'Distance ($\AA$)')
            Text(0, 0.5, r'Distance ($\AA$)')
            plt.savefig(f'{molecule_name}_all_atom_distance_array.png')
        return

    def contact_analysis(universe1, universe2, selection1, selection2, radius, molecule_name, **kwargs):
        #kwargs: beta, lambda
        #calculate all distances
        universe1 = universe1.select_atoms(selection1)
        universe2 = universe2.select_atoms(selection2)
        cont = contacts.Contacts(universe1, select=(selection1, selection2), refgroup=(universe1, universe2), radius=radius, method=method).run()
        cont_df = pd.DataFrame(cont.timeseries, columns=['Frame', 'Contacts from first Frame'])
        init_contact_matrix = cont.initial_contacts[0]
        cont_dist_ = contacts.Contacts(universe1, select=(selection1, selection2), refgroup=(universe1, universe2), radius=radius, method=method).run()
        cont_dist = pd.DataFrame(cont.timeseries, columns=['Frame', 'Contacts from first Frame'])
        reference_universe = universe1
        reference_universe.trajectory[-1]
        universe1_ = universe1.select_atoms(selection1)
        universe2_ = universe2.select_atoms(selection2)
        if kwargs.get('beta') is None:
            beta = 5.0
        if kwargs.get('lambda_') is None:
            lambda_=1.5
        cont_soft_ = contacts.Contacts(universe1, select=(selection1, selection2), refgroup=[(universe1, universe2),(universe1_, universe2_)], radius=radius, method='soft_cut', kwargs={'beta': beta, 'lambda_constant': lambda_}).run()
        cont_soft = pd.DataFrame(cont_soft_.timeseries, columns=['Frame', 'Contacts from first Frame', 'Contacts from last Frame'])
        q1q2_ = contacts.q1q2(universe1, 'name CA', radius=radius*2).run()
        q1q2 = pd.DataFrame(q1q2_.timeseries, columns=['Frame', 'Q1', 'Q2'])
        #left upper plot
        plt.subplot(2,2,1)
        cont_df.plot('Frame')
        plt.xlabel('Frame')
        plt.ylabel('Fraction of Contacts')
        plt.legend()
        #right upper plot
        plt.subplot(2,2,2)
        cont_dist.plot(x='Frame')
        plt.xlabel('Frame')
        plt.ylabel('Fraction of Contacts w/ Radius cutoff')
        plt.legend()
        #left lower plot
        plt.subplot(2,2,3)
        cont_soft.plot(x='Frame')
        plt.xlabel('Frame')
        plt.ylabel('Fraction of Contacts with soft cutoff')
        plt.legend()
        #right lower plot
        plt.subplot(2,2,4)
        q1q2.plot(x='Frame')
        plt.xlabel('Frame')
        plt.ylabel('Q1 and Q2 contacts')
        plt.legend()
        #save figure
        plt.savefig(f'{molecule_name}_contacts_all.png')
        return cont_df, cont_dist, cont_soft, q1q2
    
    def salt_bridges_in_transition(universe, molecule_name, **kwargs):
        if kwargs.get('selection_acidic') is None:
            selection_acidic = "(resname ARG LYS) and (name NH* NZ)"
            print('Assuming acidic residues like (resname ARG LYS) and (name NH* NZ)')
        if kwargs.get('selection_basic') is None:
            selection_basic = "(resname ASP GLU) and (name OE* OD*)"
            print('Assuming basic residues like (resname ASP GLU) and (name OE* OD*)')
        acidic = universe.select_atoms(selection_acidic)
        basic = universe.select_atoms(selection_basic)
        radius, timeseries = 4.5, []
        for ts in universe.trajectory:
            dist = contacts.distance_array(acidic.positions, basic.positions)
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append([ts.frame, n_contacts])
        timeseries_ = np.array(timeseries)
        timeseries_hbond = pd.DataFrame(timeseries_, Columns=['Frame', 'Number of Contacts'])
        timeseries_hbond.plot(x='Frame')
        plt.ylabel('Number of salt bridges')
        plt.savefig(f'{molecule_name}_salt_bridges.png')
        return timeseries_hbond
    
    def path_similarity(universe1, universe2, universe3, universe4, universe5, reference_universe, molecule_name, **kwargs):
        #method = Hausdorff or discrete Fréchet
        #kwargs: reference_selection, labels
        if kwargs.get('reference_selection') is None:
            reference_selection = 'name CA'
        if kwargs.get('labels') is None:
            labels = ['DCD', 'DCD2', 'XTC', 'NAMD', 'mixed']
            print:(labels)
        ps =  psa.PSAnalysis([universe1, universe2, universe3, universe4, universe5], labels=labels, reference=reference_universe, select=reference_selection, path_select='name_CA')
        ps.generate_paths(align=True, save=False, weights='mass')
        ps_ = ps
        ps.run(metric='hausdorff')
        hausdorff = ps.D
        ps_.run(metric='discrete_frechet')
        frechet = ps_.D
        #left plot (hausdorff)
        plt.subplot(2,1,1)
        hausdorff.plot_annotated_heatmap(linkage='single')
        plt.title('Hausdorff Distance/RMSD between conformations P,Q')
        plt.xlabel('P')
        plt.ylabel('Q')
        #right plot (frechet)
        plt.subplot(2,1,2)
        frechet.plot_annotated_heatmap(linkage='single')
        plt.title('Discrete Frechet/lowest possible RMSD between conformations P,Q')
        plt.xlabel('P')
        plt.ylabel('Q')
        plt.savefig(f'{molecule_name}_path_similarity.png')
        return

    def ensemble_familiarity(universe1, universe2, universe3, method, molecule_name, **kwargs):
        #methods: harmonic, clustering,  dimension reduction ens.fam.,
        #kwargs: labels, submethod, n_cores, number_of_clusters1, number_of_clusters2, num_samples1, num_samples2, eps1, eps2, universe4
        if kwargs.get('labels') is None:
             labels = ['DCD', 'DCD2', 'XTC', 'NAMD']
             print(labels)
        if method == 'harmonic':
            if kwargs.get('universe4') is None:
                raise TypeError('UniverseError: Only three universes have been given but 4 are necessary. Please specifiy universe4 as a kwarg.')
                sys.exit()
            hes, details = encore.hes([universe1, universe2, universe3, universe4], select='backbone', align=True, cov_estimator='shrinkage', weights='mass')
            fig, ax = plt.subplots()
            im = plt.imshow(hes)
            plt.xticks(np.arange(4), labels)
            plt.yticks(np.arange(4), labels)
            plt.title('Harmonic ensemble similarity')
            cbar = fig.colorbar(im)
            plt.savefig(f'{molecule_name}_harmonic_ensemble.png')
        if method == 'clustering':
            if kwargs.get('submethod') == 'default':
                labels = ['DCD', 'DCD2', 'NAMD']
                ces0, details0 = encore.ces([universe1, universe2, universe3])
                cluster_collection = details0['clustering'][0]
                first_cluster = cluster_collection.clusters[0]
                fig0, ax0 = plt.subplots()
                im0 = plt.imshow(ces0, vmax=np.log(2), vmin=0)
                plt.xticks(np.arange(3), labels)
                plt.yticks(np.arange(3), labels)
                plt.title('Clustering ensemble similarity')
                cbar0 = fig0.colorbar(im0)
                cbar0.set_label('Jensen-Shannon divergence')
            if kwargs.get('submethod') == 'affinity_propagation':
                clustering_method = clm.AffinityPropagationNative(preference=-1.0, damping=0.9, max_iter=200, convergence_iter=30, add_noise=True)
                if kwargs.get('n_cores') is None:
                    n_cores = 4
                ces1, details1 = encore.ces([universe1, universe2, universe3], select='name CA', clustering_method=clustering_method, ncores=n_cores)
                fig1, ax1 = plt.subplots()
                im1 = plt.imshow(ces1, vmax=np.log(2), vmin=0)
                plt.xticks(np.arange(3), labels)
                plt.yticks(np.arange(3), labels)
                plt.title('Clustering ensemble similarity')
                cbar1 = fig1.colorbar(im1)
                cbar1.set_label('Jensen-Shannon divergence')
            if kwargs.get('submethod') == 'k-means':
                if kwargs.get('number_of_clusters1') is None:
                    number_of_clusters1 = 6
                if kwargs.get('number_of_clusters2') is None:
                    number_of_clusters2 = 12
                if kwargs.get('min_samples1') is None:
                    min_samples1 = 5
                if kwargs.get('min_samples2') is None:
                    min_samples2 = 5
                if kwargs.get('eps1') is None:
                    eps1 = 0.5
                if kwargs.get('eps2') is None:
                    eps2 = 1.0
                km1 = clm.KMeans(number_of_clusters1, init = 'k-means++', algorithm="auto")
                km2 = clm.KMeans(number_of_clusters2, init = 'k-means++', algorithm="auto")
                db1 = clm.DBSCAN(eps=eps1, min_samples=min_samples1, algorithm='auto', leaf_size=30)
                db2 = clm.DBSCAN(eps=eps2, min_samples=min_samples2, algorithm='auto', leaf_size=30)
                ces2, details2 = encore.ces([universe1, universe2, universe3], select='name CA', clustering_method=[km1, km2, db1, db2], ncores=4)
                titles = [f'k-means w/ {number_of_clusters1} clusters', f'k-means w/ {number_of_clusters2} clusters', f'DBSCAN eps={eps1}', f'DBSCAN eps={eps2}']
                fig2, axes = plt.subplots(1, 4, sharey=True, figsize=(15, 3))
                for i, (data, title) in enumerate(zip(ces2, titles)):
                    imi = axes[i].imshow(data, vmax=np.log(2), vmin=0)
                    axes[i].set_xticks(np.arange(3))
                    axes[i].set_xticklabels(labels)
                    axes[i].set_title(title)
                plt.yticks(np.arange(3), labels)
                cbar2 = fig2.colorbar(imi, ax=axes.ravel().tolist())
                cbar2.set_label('Jensen-Shannon divergence')
            if kwargs.get('submethod') is None:
                raise TypeError('SubmethodError: If given the method clustering, this function needs a submethod. possible submethods are: \ndefault, affinity_propagation, k-means.')
            avgs, stds = encore.ces([universe1, universe2, universe3], select='name CA', clustering_method=clustering_method, estimate_error=True, ncores=4)
            np.savetxt(f'{molecule_name}_clustering_avgs.txt', avgs)
            np.savetxt(f'{molecule_name}_clustering_stds.txt', stds)
        if method == 'dimension_reduction':
            if kwargs.get('labels') is None:
                labels=['DCD', 'DCD2', 'NAMD']
            if kwargs.get('submethod') == 'default':
                dres3, details0 = encore.dres([universe1, universe2, universe3])
                reduced = details0['reduced_coordinates'][0]
                fig3, ax3 = plt.subplots()
                im3 = plt.imshow(dres3, vmax=np.log(2), vmin=0)
                plt.xticks(np.arange(3), labels)
                plt.yticks(np.arange(3), labels)
                plt.title('Dimension reduction ensemble similarity')
                cbar3 = fig3.colorbar(im0)
                cbar3.set_label('Jensen-Shannon divergence')
                plt.savefig(f'{molecule_name}_dim_red_similarity.png')
                rdfig0 = plt.figure()
                rdax0 = rdfig0.add_subplot(111, projection='3d')
                for data, label in helpers.zip_data_with_labels(reduced=reduced):
                    rdax0.scatter(*data, label=label)
                plt.legend()
                plt.savefig(f'{molecule_name}_reduced_dimensions_plot.png')
            if kwargs.get('submethod') == 'stochastic_proximity_embedding':
                dim_red_method = drm.StochasticProximityEmbeddingNative(dimension=3, min_lam=0.2, max_lam=1.0, ncycle=50, nstep=1000)
                dres1, details1 = encore.dres([universe1, universe2, universe3], select='name CA', dimensionality_reduction_method=dim_red_method, nsamples=1000, ncores=4)
                fig1, ax1 = plt.subplots()
                im1 = plt.imshow(dres1, vmax=np.log(2), vmin=0)
                plt.xticks(np.arange(3), labels)
                plt.yticks(np.arange(3), labels)
                plt.title('Dimension reduction ensemble similarity')
                cbar1 = fig1.colorbar(im1)
                cbar1.set_label('Jensen-Shannon divergence')
                plt.savefig(f'{molecule_name}_SPE_similarity.png')
            if kwargs.get('submethod') == 'PCA':
                pc1 = drm.PrincipalComponentAnalysis(dimension=1, svd_solver='auto')
                pc2 = drm.PrincipalComponentAnalysis(dimension=2, svd_solver='auto')
                pc3 = drm.PrincipalComponentAnalysis(dimension=3, svd_solver='auto')
                pc4 = drm.PrincipalComponentAnalysis(dimension=4, svd_solver='auto')
                dres2, details2 = encore.dres([universe1, universe2, universe3], select='name CA', dimensionality_reduction_method=[pc1, pc2, pc3, pc4], ncores=4)
                #plot PCA submethod results as colorboxes
                titles = ['Dim = {}'.format(n) for n in range(1, 5)]
                fig2, axes = plt.subplots(1, 4, sharey=True, figsize=(15, 3))
                for i, (data, title) in enumerate(zip(dres2, titles)):
                    imi = axes[i].imshow(data, vmax=np.log(2), vmin=0)
                    axes[i].set_xticks(np.arange(3))
                    axes[i].set_xticklabels(labels)
                    axes[i].set_title(title)
                plt.yticks(np.arange(3), labels)
                cbar2 = fig2.colorbar(imi, ax=axes.ravel().tolist())
                cbar2.set_label('Jensen-Shannon divergence')
                plt.savefig(f'{molecule_name}_PCA_similarity.png')
                #plot stochastic proximity embedding as violin
                rd_p1, rd_p2, rd_p3, _ = details2['reduced_coordinates']
                rd_p1_fig, rd_p1_ax = plt.subplots(figsize=(4, 8))
                split_data = [x[0].reshape((-1,)) for x in zip_data_with_labels(rd_p1)]
                rd_p1_ax.violinplot(split_data, showextrema=False)
                rd_p1_ax.set_xticks(np.arange(1, 4))
                rd_p1_ax.set_xticklabels(labels)
                plt.savefig(f'{molecule_name}_PCA_violin.png')
                #plot stochastic proximity embedding in 3D
                rd_p3_fig = plt.figure(figsize=(8, 6))
                rd_p3_ax = rd_p3_fig.add_subplot(111, projection='3d')
                for data, label in zip_data_with_labels(rd_p3):
                    rd_p3_ax.scatter(*data, label=label)
                rd_p3_ax.set_xlabel('PC 1')
                rd_p3_ax.set_ylabel('PC 2')
                rd_p3_ax.set_zlabel('PC 3')
                plt.legend()
                plt.savefig(f'{molecule_name}_PCA_3D.png')
                # save errors
                avgs, stds = encore.dres([universe1, universe2, universe3], select='name CA', dimensionality_reduction_method=dim_red_method, estimate_error=True, ncores=4)
                np.savetxt(f'{molecule_name}_clustering_avgs.txt', avgs)
                np.savetxt(f'{molecule_name}_clustering_stds.txt', stds)
        return 
    
    def eval_convergence(universe, molecule_name, **kwargs):
        #kwargs: method:{'default', 'k-means', 'dim-reduction', 'PCA', 'all'}
        #number_of_clusters1, number_of_clusters2, number_of_clusters3, dim1, dim2, dim3
        if kwargs.get('window_size') is None:
            window_size = 10
        if kwargs.get('method') == 'default': 
            ces_conv = encore.ces_convergence(universe, window_size, select='name CA')
            ces_fig, ces_ax = plt.subplots()
            plt.plot(ces_conv)
            ces_ax.set_xlabel('Window')
            ces_ax.set_ylabel('Jensen-Shannon divergence')
            plt.savefig(f'{molecule_name}_JS_divergence.png')
        if kwargs.get('method') == 'k-means':
            if kwargs.get('number_of_clusters1') is None:
                raise('K-MeansError: please specify an amount of clusters to be used for k-means analysis. Necessary for this are:\nnumber_of_clusters1, number_of_clusters2, number_of_clusters3')
            km1 = clm.KMeans(number_of_clusters1, init = 'k-means++', algorithm="auto")
            km2 = clm.KMeans(number_of_clusters2, init = 'k-means++', algorithm="auto")
            km3 = clm.KMeans(number_of_clusters3, init = 'k-means++', algorithm="auto")
            ces_conv = encore.ces_convergence(universe, window_size, select='name CA', clustering_method=[km1, km2, km3])
            labels = [f'{number_of_clusters1} clusters', f'{number_of_clusters2} clusters', f'{number_of_clusters3} clusters']
            ces_fig2, ces_ax2 = plt.subplots()
            for data, label in zip(ces_conv2.T, labels):
                plt.plot(data, label=label)
            ces_ax2.set_xlabel('Window')
            ces_ax2.set_ylabel('Jensen-Shannon divergence')
            plt.legend()
            plt.savefig(f'{molecule_name}_kMeans_convergence.png')
        if kwargs.get('method') == 'dim-reduction':
            dim_red = encore.dres_convergence(universe, 10, select='name CA')
            dres_fig, dres_ax = plt.subplots()
            plt.plot(dim_red)
            dres_ax.set_xlabel('Window')
            dres_ax.set_ylabel('Jensen-Shannon divergence')
            plt.savefig(f'{molecule_name}_dim_red_ens_sim.png')
        if kwargs.get('method') == 'PCA':
            pc1 = drm.PrincipalComponentAnalysis(dimension=dim1, svd_solver='auto')
            pc2 = drm.PrincipalComponentAnalysis(dimension=dim2, svd_solver='auto')
            pc3 = drm.PrincipalComponentAnalysis(dimension=dim3, svd_solver='auto')
            dres_conv2 = encore.dres_convergence(universe, 10, select='name CA', dimensionality_reduction_method=[pc1, pc2, pc3])
            labels=[f'{dim1}D', f'{dim2}D', f'{dim3}D']
            dres_fig2, dres_ax2 = plt.subplots()
            for data, label in zip(dres_conv2.T, labels):
                plt.plot(data, label=label)
            dres_ax2.set_xlabel('Window')
            dres_ax2.set_ylabel('Jensen-Shannon divergence')
            plt.legend()
            plt.savefig(f'{molecule_name}_pca_convergence.png')
        if kwargs.get('method') == 'all':
            ces_conv = encore.ces_convergence(universe, window_size, select='name CA')
            km1 = clm.KMeans(number_of_clusters1, init = 'k-means++', algorithm="auto")
            km2 = clm.KMeans(number_of_clusters2, init = 'k-means++', algorithm="auto")
            km3 = clm.KMeans(number_of_clusters3, init = 'k-means++', algorithm="auto")
            ces_conv = encore.ces_convergence(universe, window_size, select='name CA', clustering_method=[km1, km2, km3])
            dim_red = encore.dres_convergence(universe, window_size, select='name CA')
            pc1 = drm.PrincipalComponentAnalysis(dimension=dim1, svd_solver='auto')
            pc2 = drm.PrincipalComponentAnalysis(dimension=dim2, svd_solver='auto')
            pc3 = drm.PrincipalComponentAnalysis(dimension=dim3, svd_solver='auto')
            dres_conv2 = encore.dres_convergence(universe, window_size, select='name CA', dimensionality_reduction_method=[pc1, pc2, pc3])
            fig, ax = plt.subplots()
            labels = [f'{number_of_clusters1} clusters', f'{number_of_clusters2} clusters', f'{number_of_clusters3} clusters']
            labels_=[f'{dim1}D', f'{dim2}D', f'{dim3}D']
            plt.plot(ces_conv, label='clustering ensemble')
            for data, label in zip(ces_conv2.T, labels):
                plt.plot(data, label=label)
            plt.plot(dim_red, label='dim-red ensemble similarity')
            for data, label in zip(dres_conv2.T, labels_):
                plt.plot(data, label=label)
            dres_ax2.set_xlabel('Window')
            dres_ax2.set_ylabel('Jensen-Shannon divergence')
            plt.legend()
            plt.savefig(f'{molecule_name}_convergence_all.png')
        return
    
    def elastic_network(universe1, universe2, molecule_name, method, **kwargs):
        if method == 'default':
            nma1 = gnm.GNMAnalysis(universe1, select='name CA', cutoff=7.0)
            nma1.run()
            nma2 = gnm.GNMAnalysis(universe2, select='name CA', cutoff=7.0)
            nma2.run()
            eigenvalues1 = [res[1] for res in nma1.results]
            eigenvalues2 = [res[1] for res in nma2.results]
            histfig, histax = plt.subplots(nrows=2, sharex=True, sharey=True)
            histax[0].hist(eigenvalues1)
            histax[1].hist(eigenvalues2)
            histax[1].set_xlabel('Eigenvalue')
            histax[0].set_ylabel('Frequency')
            histax[1].set_ylabel('Frequency')
            plt.savefig(f'{molecule_name}_elastic_net_eigenvalues.png')
            time1 = [res[0] for res in nma1.results]
            time2 = [res[0] for res in nma2.results]
            linefig, lineax = plt.subplots()
            plt.plot(time1, eigenvalues1, label='DCD')
            plt.plot(time2, eigenvalues2, label='DCD2')
            lineax.set_xlabel('Time (ps)')
            lineax.set_ylabel('Eigenvalue')
            plt.legend()
            plt.savefig(f'{molecule_name}_elastic_net_eigen_variation.png')
        if method == 'close_contact':
            nma_close = gnm.closeContactGNMAnalysis(universe1, select='name CA', cutoff=7.0, weights='size')
            nma_close.run()
            eigenvalues_close = [res[1] for res in nma_close.results]
            plt.hist(eigenvalues_close)
            plt.xlabel('Eigenvalue')
            plt.ylabel('Frequency')
            plt.savefig(f'{molecule_name}_elastic_net_close_contacts.png')
            time_close = [res[0] for res in nma_close.results]
            ax = plt.plot(time_close, eigenvalues_close)
            plt.xlabel('Time (ps)')
            plt.ylabel('Eigenvalue')
            plt.savefig(f'{molecule_name}_elastic_net_close_variation.png')
        else: 
            raise TypeError('method needs to be either default or close_contact.')
        return
    
    def average_rdf(universe, molecule_name, **kwargs):
        #kwargs: exclusion=(True/False), site_specific=(density/no_density), selection1, selection2, selection3
        #method=(atom_to_atom: atom_group_pair, atom_from, atom_to)
        if kwargs.get('exclusion') == False:
            sele = universe.select_atoms(selection)
            waterverse = universe.select_atoms('resname SOL')
            irdf = rdf.InterRDF(sele, waterverse, nbins=75,  range=(0.0, 15.0))
            irdf.run()
            plt.plot(irdf.bins, irdf.rdf)
            plt.xlabel('Radius (angstrom)')
            plt.ylabel('Radial distribution')
            plt.savefig(f'{molecule_name}_radial_distribution_NoExclusion.png')
        if kwargs.get('exclusion') == True:
            irdf2 = rdf.InterRDF(sele, sele, exclusion_block=(1, 1))
            irdf2.run()
            plt.plot(irdf2.bins, irdf2.rdf)
            plt.xlabel('Radius (angstrom)')
            plt.ylabel('Radial distribution')
            plt.savefig(f'{molecule_name}_radial_distribution_Exclusion.png')
        if kwargs.get('site_specific') == 'density':
            if kwargs.get('selection1') is None:
                selection1 = 'resid50'
                print('No selection was given, calculating with Residue 50 as standard value.')
            sele = universe.select_atoms(selection1)
            waterverse = universe.select_atoms('resname SOL')
            irdf = rdf.InterRDF(sele, water, nbins=75,  range=(0.0, 15.0))
            irdf.run()
            plt.plot(irdf.bins, irdf.rdf)
            plt.xlabel('Radius (angstrom)')
            plt.ylabel('Radial distribution')
            plt.savefig(f'{molecule_name}_rdf_density.png')
        if kwargs.get('site_specific') == 'no_density':
            if kwargs.get('selection1') is None:
                selection1 = 'resname THR'
            if kwargs.get('selection2') is None:
                selection2 = selection1
                print('No selection 2 given, assuming selection2=selection1, and calculating the average RDF for selection1-atom group self-overlap.')
            irdf2 = rdf.InterRDF(selection1, selection2, exclusion_block=(1, 1))
            irdf2.run()
            plt.plot(irdf2.bins, irdf2.rdf)
            plt.xlabel('Radius (angstrom)')
            plt.ylabel('Radial distribution')
            plt.savefig((f'{molecule_name}_rdf_NoDensity.png'))
        if kwargs.get('method') == 'atom_to_atom':
            if kwargs.get('selection1') is None:
                selection1 = 'resid 60 and name CA'
            if kwargs.get('selection2') is None:
                selection2 = 'resid 61 and name CA'
            if kwargs.get('selection3') is None:
                selection3 = 'resid 62 and name CA'
            atom1 = universe.select_atoms(selection1)
            atom2 = universe.select_atoms(selection2)
            atom3 = universe.select_atoms(selection3)
            waterverse = universe.select_atoms('resname SOL and sphzone 15 group sel_a', sel_a=ca60)
            ags = [[ca60+ca61, water], [ca62, water]]
            ss_rdf = rdf.InterRDF_s(universe, ags, nbins=75, range=(0.0, 15.0), density=True)
            ss_rdf.run()
            if kwargs.get('atom_group_pair') is None:
                atom_group_pair = 0
            if kwargs.get('atom_from') is None:
                atom_from = 1
                print('No atom to calculate the distance from has been given, assuming you want to calculate the distance from atom 1.')
            if kwargs.get('atom_to') is None:
                atom_to = 2
                print('No atom to calculate the distance to has been given, assuming you want to calculate the distance towards atom 2.')
            atom1_atom2_distance = ss_rdf.rdf[atom_group_pair][atom_from][atom_to]
            plt.xlabel('Radius (angstrom)')
            plt.ylabel('Radial distribution')
            plt.title(f'RDF between atoms {atom_from} and {atom_to} in group {atom_group_pair}')
            plt.savefig(f'{molecule_name}_dist_atom{atom_from}_atom{atom_to}.png')
        return 
    
    def dihedral_analysis(universe, molecule_name, **kwargs):
        #kwargs: verbose, angle, upper, lower, show_reference_angles
        if kwargs.get('show_reference_angles') is None:
            show_reference_angles = True
        if kwargs.get('lower') is None:
            lower = np.random.randint(0,len(protein.residues))
        if kwargs.get('upper') is None:
            upper = lower + 2
        protein = universe.select_atoms('protein')
        #calculate angles
        for res in universe.residues[:]:
            if kwargs.get('angle') is None:
                angle = res.phi_selection()
            if kwargs.get('angle') == 'phi':
                angle = res.phi_selection()
            if kwargs.get('angle') == 'psi':
                angle = res.psi_selection()
            if kwargs.get('angle') == 'omega':
                angle = res.omega_selection()
            if kwargs.get('angle') == 'chi1':
                angle = res.chi1_selection()
            if angles is None:
                names = None
            else:
                names = angles.names
            if kwargs.get('verbose') is not None:
                print(f' Currently printed angle is {angle}.')
                print(f'{res.resname}: {names} ')
        #calculate special angles (omegas, then dihedrals from that)
        omegas = [res.omega_selection() for res in protein.residues[lower:upper]]
        dihedral = dihedrals.Dihedral(omegas).run()
        #calculate Ramachandran
        rama = dihedrals.Ramachandran(protein).run()
        #calculate Janin
        janin = dihedrals.Janin(protein).run()
        #plot dihedrals
        labels = [f'Res {n}' for n in np.arange(lower, upper)]
        for ang, label in zip(dihedral.angles.T, labels):
            plt.plot(ang, label=label)
        plt.xlabel('Frame')
        plt.ylabel('Angle [°]')
        plt.legend()
        plt.savefig(f'{molecule_name}_dihedral_angles_res{lower}_res{upper}.png')
        #plot polar plot of dihedrals
        fig_polar = plt.figure()
        ax_polar = fig_polar.add_subplot(111, projection='polar')
        frames = np.arange(10)
        for res, label in zip(dihedral.angles.T, labels):
            c = ax_polar.plot(res, frames, label=label)
        plt.legend()
        plt.savefig(f'{molecule_name}_polar_dihedral_angles_res{lower}_res{upper}.png')
        #Ramachandran Plot
        rama.plot(color='black', marker='.', ref=show_reference_angles)
        plt.savefig(f'{molecule_name}_ramachandran.png')
        #Janin Plot
        janin.plot(ref=show_reference_angles, marker='.', color='black')
        plt.savefig(f'{molecule_name}_janin.png')
        return

    def MD_principal_component_analysis(molecule_name, universe, **kwargs):
        aligner = align.AlignTraj(universe, universe, select="backbone", in_memory=True).run()
        pc = pca.PCA(universe, select="backbone", align=True, mean=None, n_components=None).run()
        backbone = u.select_atoms("backbone")
        n_bb = len(backbone)
        #pc.variance[0] gives the variance of the first variable
        if kwargs.get('lower') is None:
            lower = 0
        if kwargs.get('upper') is None:
            upper = n_bb
        plt.plot(pc.cumulated_variance[lower:upper])
        plt.xlabel("Principal component")
        plt.ylabel('Cumulative Variance')
        transformed = pc.transform(backbone, n_components=5)
        df = pd.DataFrame(transformed, columns=[f'PC{i+1}' for i in range(5)])
        df["Time (ps)"] = df.index * u.trajectory.dtn
        g = sns.PairGrid(df, hue="Time (ps)", palette=sns.color_palette("Oranges_d", n_colors=len(df)))
        g.map(plt.scatter, marker=".")
        plt.savefig(f'{molecule_name}_pca_transformed.png')
        pc1 = pc.p_components[:, 0]
        trans1 = transformed[:, 0]
        projected = np.outer(trans1, pc1) + pc.mean
        coordinates = projected.reshape(len(trans1), -1, 3)
        proj1 = mda.Merge(backbone)
        proj1.load_new(coordinates, order="fac")
        movie = MovieMaker(view, output='pc1.gif', in_memory=True)
        movie.make()
        for i in range(5):
            cc = pca.cosine_content(transformed, i)
            print(f'Cosine content for PC {i+1} = {cc:.3f}')
        melted = pd.melt(df, id_vars=["Time (ps)"], var_name="PC", value_name="Value")
        g = sns.FacetGrid(melted, col="PC")
        g.map(sns.lineplot, "Time (ps)", "Value", ci=None)
        plt.savefig(f'{molecule_name}_pca_multi.png')
        return
    
    def diffusion_map(molecule_name, universe, **kwargs):
        dmap = diffusionmap.DiffusionMap(u, select='backbone', epsilon=2)
        dmap.run()
        fig, ax = plt.subplots()
        ax.plot(dmap.eigenvalues[1:16])
        plt.savefig(f'{molecule_name}_diffusion_map_eigenv.png')
        transformed = dmap.transform(5, time=1)
        df = pd.DataFrame(transformed, columns=['Mode{}'.format(i+2) for i in range(5)]) 
        df['Time (ps)'] = df.index * u.trajectory.dt
        g = sns.PairGrid(df, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df)))
        g.map(plt.scatter, marker='.')
        plt.savefig(f'{molecule_name }_diffusion_map_scatter.png')
        return

    def persistence_length(molecule_name, universe, **kwargs):
        chains = universe.atoms.fragments
        backbones = [ch.select_atoms('not name O* H*') for ch in chains]
        sorted_backbones = [polymer.sort_backbone(bb) for bb in backbones]
        persistence_length = polymer.PersistenceLength(sorted_backbones)
        persistence_length.run()
        plen.plot()
        plt.savefig(f'{molecule_name}_persistence_length.png')
        return
    
    def HOLE(molecule_name, universe, universe_dim, **kwargs):
        #Default Random Seeds von menschlichem Random Seed Generator Anna T.
        if universe_dim == 'pdb':
            if kwargs.get('rand') is None:
                rand = 72843
            profiles = hole2.hole(universe, executable='~/hole2/exe/hole', outfile='hole.out', sphpdb_file='hole.sph', vdwradii_file=None, random_seed=rand)
            rxn_coords = profiles[0].rxn_coord 
            hole2.create_vmd_surface(filename='hole.vmd', sphpdb='hole.sph', sph_process='~/hole2/exe/sph_process')
        if universe_dim == 'traj':
            #Pore Radius
            ha = hole2.HoleAnalysis(universe, select='protein', cpoint='center_of_geometry', executable='~/hole2/exe/hole')
            if kwargs.get('rand') is None:
                rand = 32456
            ha.run(random_seed=rand)
            gathered = ha.gather()
            flat = ha.gather(flat=True)
            if kwargs.get('bins') is None:
                bins=100
            radii, edges = ha.bin_radii(bins=bins, range=None)
            means, edges_ = ha.histogram_radii(bins=bins, range=None, aggregator=np.mean)
            midpoints=0.5*(edges[1:]+edges[:-1])
            plt.plot(midpoints, means)
            plt.ylabel(r"Mean HOLE radius $R$ ($\AA$)")
            plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
            plt.savefig(f'{molecule_name}_pore_coordinate.png')
            #Minimum Pore Radius
            min_radii = ha.min_radius()
            plt.plot(min_radii[:,0], min_radii[:,1])
            plt.ylabel('Minimum HOLE radius $R$ ($\AA$)')
            plt.xlabel('Frame')
            plt.savefig(f'{molecule_name}_minimum_pore_radius.png')
            #Create VMD Video of pore over time
            ha.create_vmd_surface(filename='holeanalysis.vmd')
            if kwargs.get('3D') == True:
                ha.plot3D()
                ha.plot_mean_profile(bins=bins, n_std=1, color='blue', fill_alpha=0.2, legend=True)
        return
    
    def mass_density(molecule_name, universe, **kwargs):
        density = lin.LinearDensity(universe.atoms, grouping='atoms').run()
        plt.plot(np.linspace(0,50,200), density.results['x']['pos'], density.results['x']['pos_std'], legend=['Position X', 'StD of Position X'])
        plt.plot(np.linspace(0,50,200), density.results['y']['pos'], density.results['y']['pos_std'], legend=['Position Y', 'StD of Position Y'])
        plt.plot(np.linspace(0,50,200), density.results['z']['pos'], density.results['z']['pos_std'], legend=['Position Z', 'StD of Position Z'])
        plt.savefig(f'{molecule_name}_mass_density.png')
        return
    
    def solvent_density(molecule_name, universe, waterverse, **kwargs):
        workflow = [trans.unwrap(u.atoms), 
            trans.center_in_box(protein, center='geometry'),
            trans.wrap(water, compound='residues'), 
            trans.fit_rot_trans(protein, protein, weights='mass')]
        universe.trajectory.add_transformations(*workflow)
        view = nv.show_mdanalysis(universe)
        view.add_representation('point', selection='resname SOL')
        view.render_image()
        OW = universe.select_atoms('name OW')
        density_ = density.DensityAnalysis(OW, delta=4.0, padding=2)
        density_.run()
        grid = density_.density.grid
        density_.density.convert_density('TIP4P')
        #from skimage import measure, from mpl_toolkits.mpl3d import Axes3D
        if kwargs.get('iso_val') is None:
            iso_val = 0.5
        if kwargs.get('method') == 'default':
            verts, faces, _, _ = measure.marching_cubes_lewiner(density_.density.grid, iso_val, spacing=density_.density.delta)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2], lw=1, alpha=0.1)
            plt.savefig(f'{molecule_name}_solvent_density.png')
        #import pyvista as pv
        if kwargs.get('method') == 'advanced':
            pv.set_plot_theme("document")
            x, y, z = np.meshgrid(mx, my, mz, indexing="ij")
            mesh = pv.StructuredGrid(x, y, z)
            mesh["density"] = dens.density.grid.T.flatten() # note transpose
            contours = mesh.contour([0.5, 1.2])
            p = pv.Plotter(notebook=True)
            p.background_color = 'white'
            p.add_mesh(mesh.outline(), color="k")  # box lines
            p.add_mesh(contours, opacity=0.2)  # surfaces
            p.savefig(f'{molecule_name}_solvent_density.png')
        return
    
    def workflow():
        return
