# coding=utf-8

import kbp2
import re
import numpy as np




if __name__ == '__main__':

    dict1= {'HSD-326_ACHA':{'NE2':'N32',\
                            'ND1':'N27',\
                            'HD1':'H28',\
                            'CE1':'C30',\
                            'CD2':'C33',\
                            'CG':'C29',\
                            'HD2':'H34',\
                            'HE1':'H31',\
                            'CB':'H99'},\

            'HSD-325_ACHA':{'NE2':'N24',\
                            'ND1':'N19',\
                            'HD1':'H20',\
                            'CE1':'C22',\
                            'CD2':'C25',\
                            'CG':'C21',\
                            'HD2':'H26',\
                            'HE1':'H23',\
                            'CB':'H101'},\

            'HSE-276_ACHA':{'NE2':'N5',\
                            'ND1':'N1',\
                            'HE1':'H4',\
                            'CE1':'C3',\
                            'CD2':'C6',\
                            'CG':'C2',\
                            'HD2':'H7',\
                            'CB':'H100'},\

            'HSD-411_ACHA':{'NE2':'N40',\
                            'ND1':'N35',\
                            'HD1':'H36',\
                            'CE1':'C38',\
                            'CD2':'C41',\
                            'CG':'C37',\
                            'HD2':'H42',\
                            'HE1':'H39',\
                            'CB':'H98'},\

            'TYR-280_ACHA':{'CE2':'C18',\
                            'CZ':'C13',\
                            'CE1':'C11',\
                            'CD1':'C9',\
                            'CG':'C8',\
                            'CD2':'C16',\
                            'OH':'O14',\
                            'HH':'H15',\
                            'HE1':'H12',\
                            'HD1':'H10',\
                            'HD2':'H17',\
                            'C3':'H102'},\

            'HEM-2_EHEM':{'HMA':'H106',\
                          'OMA':'O72',\
                          'CMA':'C71',\
                          'C3A':'C50',\
                          'C4A':'C51',\
                          'C2A':'C49',\
                          'C1A':'C48',\
                          'NA':'N44',\
                          'CHB':'C65',\
                          'HB':'H66',\
                          'C1B':'C52',\
                          'NB':'N45',\
                          'C2B':'C53',\
                          'CMB':'C73',\
                          'HMB1':'H74',\
                          'HMB2':'H75',\
                          'HMB3':'H76',\
                          'C3B':'C54',\
                          'C11':'C90',\
                          'H1A':'H91',\
                          'O11':'O92',\
                          'HO11':'H93',\
                          'C12':'H97',\
                          'C4B':'C55',\
                          'CHC':'C67',\
                          'HC':'H68',\
                          'C1C':'C56',\
                          'C2C':'C57',\
                          'CMC':'C77',\
                          'HMC1':'H78',\
                          'HMC2':'H79',\
                          'HMC3':'H80',\
                          'C3C':'C58',\
                          'CAC':'C81',\
                          'HAC':'H82',\
                          'CBC':'C83',\
                          'HBC1':'H84',\
                          'HBC2':'H85',\
                          'C4C':'C59',\
                          'NC':'N46',\
                          'CHD':'C69',\
                          'HD':'H70', \
                          'C1D': 'C60', \
                          'C2D': 'C61', \
                          'CMD': 'C86', \
                          'HMD1': 'H87', \
                          'HMD2': 'H88', \
                          'HMD3': 'H89', \
                          'C3D': 'C62', \
                          'C4D': 'C63', \
                          'ND': 'N47', \
                          'CHA':'C64',\
                          'HA':'H105',\
                          'FE':'Fe43'},\

            'CU1-1_META':{'CU':'Cu94'},\

            'PER-1_PERI':{'OH1':'O95',\
                          'OH2':'O96'}
            }

    # for residue_descr, atoms in dict1.iteritems():
    #     resname, resid, segname = re.split('[-_]', residue_descr)
    #     for atom in atoms.keys():
    #         print '.or. (resid %s .and. segid %s .and. type %s) -' % (resid,segname,atom)

#     string = """define small sele ((resid 325 .and. segid ACHA .and. type HE1) -
# .or. (resid 325 .and. segid ACHA .and. type CD2) -
# .or. (resid 325 .and. segid ACHA .and. type CE1) -
# .or. (resid 325 .and. segid ACHA .and. type HD1) -
# .or. (resid 325 .and. segid ACHA .and. type ND1) -
# .or. (resid 325 .and. segid ACHA .and. type CG) -
# .or. (resid 325 .and. segid ACHA .and. type NE2) -
# .or. (resid 325 .and. segid ACHA .and. type HD2) -
# .or. (resid 280 .and. segid ACHA .and. type CE1) -
# .or. (resid 280 .and. segid ACHA .and. type CG) -
# .or. (resid 280 .and. segid ACHA .and. type CZ) -
# .or. (resid 280 .and. segid ACHA .and. type HH) -
# .or. (resid 280 .and. segid ACHA .and. type CE2) -
# .or. (resid 280 .and. segid ACHA .and. type CD2) -
# .or. (resid 280 .and. segid ACHA .and. type HE1) -
# .or. (resid 280 .and. segid ACHA .and. type HD2) -
# .or. (resid 280 .and. segid ACHA .and. type HD1) -
# .or. (resid 280 .and. segid ACHA .and. type OH) -
# .or. (resid 280 .and. segid ACHA .and. type C3) -
# .or. (resid 280 .and. segid ACHA .and. type CD1) -
# .or. (resid 411 .and. segid ACHA .and. type HE1) -
# .or. (resid 411 .and. segid ACHA .and. type CD2) -
# .or. (resid 411 .and. segid ACHA .and. type CE1) -
# .or. (resid 411 .and. segid ACHA .and. type HD1) -
# .or. (resid 411 .and. segid ACHA .and. type ND1) -
# .or. (resid 411 .and. segid ACHA .and. type CG) -
# .or. (resid 411 .and. segid ACHA .and. type NE2) -
# .or. (resid 411 .and. segid ACHA .and. type HD2) -
# .or. (resid 2 .and. segid EHEM .and. type HMC1) -
# .or. (resid 2 .and. segid EHEM .and. type HMC2) -
# .or. (resid 2 .and. segid EHEM .and. type HMC3) -
# .or. (resid 2 .and. segid EHEM .and. type FE) -
# .or. (resid 2 .and. segid EHEM .and. type HO11) -
# .or. (resid 2 .and. segid EHEM .and. type NB) -
# .or. (resid 2 .and. segid EHEM .and. type NC) -
# .or. (resid 2 .and. segid EHEM .and. type ND) -
# .or. (resid 2 .and. segid EHEM .and. type C2B) -
# .or. (resid 2 .and. segid EHEM .and. type C2C) -
# .or. (resid 2 .and. segid EHEM .and. type C2A) -
# .or. (resid 2 .and. segid EHEM .and. type C2D) -
# .or. (resid 2 .and. segid EHEM .and. type HMB3) -
# .or. (resid 2 .and. segid EHEM .and. type HMB2) -
# .or. (resid 2 .and. segid EHEM .and. type CHB) -
# .or. (resid 2 .and. segid EHEM .and. type CHC) -
# .or. (resid 2 .and. segid EHEM .and. type CHD) -
# .or. (resid 2 .and. segid EHEM .and. type HAC) -
# .or. (resid 2 .and. segid EHEM .and. type CMB) -
# .or. (resid 2 .and. segid EHEM .and. type HBC2) -
# .or. (resid 2 .and. segid EHEM .and. type HBC1) -
# .or. (resid 2 .and. segid EHEM .and. type OMA) -
# .or. (resid 2 .and. segid EHEM .and. type O11) -
# .or. (resid 2 .and. segid EHEM .and. type C3D) -
# .or. (resid 2 .and. segid EHEM .and. type C3A) -
# .or. (resid 2 .and. segid EHEM .and. type C3C) -
# .or. (resid 2 .and. segid EHEM .and. type C3B) -
# .or. (resid 2 .and. segid EHEM .and. type CHA) -
# .or. (resid 2 .and. segid EHEM .and. type CMD) -
# .or. (resid 2 .and. segid EHEM .and. type CMC) -
# .or. (resid 2 .and. segid EHEM .and. type HMB1) -
# .or. (resid 2 .and. segid EHEM .and. type CMA) -
# .or. (resid 2 .and. segid EHEM .and. type HB) -
# .or. (resid 2 .and. segid EHEM .and. type C12) -
# .or. (resid 2 .and. segid EHEM .and. type CAC) -
# .or. (resid 2 .and. segid EHEM .and. type HC) -
# .or. (resid 2 .and. segid EHEM .and. type C11) -
# .or. (resid 2 .and. segid EHEM .and. type HA) -
# .or. (resid 2 .and. segid EHEM .and. type HD) -
# .or. (resid 2 .and. segid EHEM .and. type C4D) -
# .or. (resid 2 .and. segid EHEM .and. type C4A) -
# .or. (resid 2 .and. segid EHEM .and. type C4B) -
# .or. (resid 2 .and. segid EHEM .and. type C4C) -
# .or. (resid 2 .and. segid EHEM .and. type H1A) -
# .or. (resid 2 .and. segid EHEM .and. type HMA) -
# .or. (resid 2 .and. segid EHEM .and. type HMD1) -
# .or. (resid 2 .and. segid EHEM .and. type HMD3) -
# .or. (resid 2 .and. segid EHEM .and. type HMD2) -
# .or. (resid 2 .and. segid EHEM .and. type CBC) -
# .or. (resid 2 .and. segid EHEM .and. type C1C) -
# .or. (resid 2 .and. segid EHEM .and. type C1B) -
# .or. (resid 2 .and. segid EHEM .and. type C1A) -
# .or. (resid 2 .and. segid EHEM .and. type C1D) -
# .or. (resid 2 .and. segid EHEM .and. type NA) -
# .or. (resid 276 .and. segid ACHA .and. type HE1) -
# .or. (resid 276 .and. segid ACHA .and. type CG) -
# .or. (resid 276 .and. segid ACHA .and. type CE1) -
# .or. (resid 276 .and. segid ACHA .and. type CD2) -
# .or. (resid 276 .and. segid ACHA .and. type ND1) -
# .or. (resid 276 .and. segid ACHA .and. type NE2) -
# .or. (resid 276 .and. segid ACHA .and. type HD2) -
# .or. (resid 326 .and. segid ACHA .and. type HE1) -
# .or. (resid 326 .and. segid ACHA .and. type CD2) -
# .or. (resid 326 .and. segid ACHA .and. type CE1) -
# .or. (resid 326 .and. segid ACHA .and. type HD1) -
# .or. (resid 326 .and. segid ACHA .and. type ND1) -
# .or. (resid 326 .and. segid ACHA .and. type CG) -
# .or. (resid 326 .and. segid ACHA .and. type NE2) -
# .or. (resid 326 .and. segid ACHA .and. type HD2) -
# .or. (resid 1 .and. segid META .and. type CU1) -
# .or. (resid 1 .and. segid PERI .and. type OH1) -
# .or. (resid 1 .and. segid PERI .and. type OH2) -
# .or. segid SWAT .or segid WAT. end
# """

    dict2= {}
    # xyz_file = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/free_step_8.xyz'
    xyz_file = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/free_step_10.xyz'
    f = open(xyz_file, 'r')
    set = 1
    i = 0
    read = False
    for line in f:
        if 'Geometry' in line:
            i += 1
            y = 1
            continue
        if (i==set) and (i !=0):
            parsed_line = line.split()
            if len(parsed_line) == 4:
                atom_name = parsed_line[0]+str(y)
                atom_coords = []
                atom_coords.append(float(parsed_line[1]))
                atom_coords.append(float(parsed_line[2]))
                atom_coords.append(float(parsed_line[3]))
                dict2[atom_name] = atom_coords
                y+=1


    pdb = kbp2.file_parser.Simple_struct_parser()
    pdb.read_pdb('/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/read_first/cco_o2_out.pdb')
    # pdb.read_crd('/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/read_first/cco_o2_out.crd')
    pdb.create_struct()
    for seg in pdb.struct.iter_segments():
            for res in seg.iter_residues():
                for residue_descr in dict1.keys():
                    resname, resid, segname = re.split('[-_]', residue_descr)
                    if (res.resname == resname) and (res.resid == int(resid)) and (res.segname == segname):
                        for atom_name in dict1[residue_descr].keys():
                            for atm in res.iter_atoms():
                                if atm['name'] == atom_name:
                                    new_coord = dict2[dict1[residue_descr][atom_name]]
                                    new_coord = np.array(new_coord)
                                    atm['coord']= new_coord
    pdb.create_struct()
    pdb_w = pdb.copy(segname=['SWAT', 'WAT'], exclude=True)
    # pdb_w.write_pdb('/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/insert_free_8.pdb')
    pdb_w.write_pdb('/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/split_new_coord_CB/cco_o2.pdb')
    # pdb_w.write_crd('/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/insert_free_8.crd')
