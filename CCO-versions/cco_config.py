# coding=utf-8

import os
import shutil
import shutil
import kbp2
import sys


def copy_patch(patch):

     patch_name = patch.keys()[0]
     patch_values = patch[patch_name]
     new_patch = {}
     new_patch[patch_name] = list(patch_values)
     return new_patch

class COX_settings():


    '''Class that contains COX information. Hardcoded patches and top/par in relation to the state and specie
    All else, like special residues should be defined later.
    It is meant to work with pka_calculation and charmm_module'''

    def __init__(self, specie = '', cycle_state = ''):

         self.specie = specie
         self.cycle_state = cycle_state

         # charge patches are always before autogenerate angle dihedral
         self.charge_patches = []

         # charge patches are always after autogenerate angle dihedral
         self.bond_patches = []

         self.top = []
         self.par =[]

         if cycle_state:
             if specie:
                self.set_information()


    def copy(self):

         new =  COX_settings()


         new.cycle_state = self.cycle_state
         new.specie = self.specie

         for top in self.top:
             new.top.append(top)

         for par in self.par:
             new.par.append(par)

         if self.charge_patches:
             for patch in self.charge_patches:
                new_patch = copy_patch(patch)
                new.charge_patches.append(new_patch)

         if self.bond_patches:
             for patch in self.bond_patches:
                new_patch = copy_patch(patch)
                new.bond_patches.append(new_patch)

         return new

    def set_information(self):

         specie = self.specie
         cycle_state = self.cycle_state


         self.top.append("/user/jdragelj/projects/cco_fluorescin/top_alw.inp")
         self.top.append("/scratch/scratch/awoelke/md_cco/toppar/patches.rtf")
         self.top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_lipid.rtf")
         self.top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_cgenff.rtf")

         self.par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all22_prot_plus_heme_and_Cu_kb.inp")
         self.par.append("/scratch/scratch/awoelke/md_cco/toppar/patches.prm")
         self.par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_lipid.prm")
         self.par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_cgenff.prm")


         ####################
         ### COMBINATIONS ###
         ####################

                                                                     ####################
                                                                     ###  PARACOCCUS  ###
                                                                     ####################

         if specie == 'A_paraccocus' and cycle_state == 'O':

            ### CHARGE PATCHES
            self.charge_patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
            self.charge_patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
            self.charge_patches.append({'CB4' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
            self.charge_patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})

            ### BOND PATCHES
            self.bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
            self.bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
            self.bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
            self.bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

         if specie == 'A_paraccocus' and cycle_state == 'R':

            ### CHARGE PATCHES
            self.charge_patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
            self.charge_patches.append({'CA11' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
            self.charge_patches.append({'CB1T' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
            self.charge_patches.append({'A3W2' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})

            ### BOND PATCHES
            self.bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
            self.bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'EISW' : ['HOH-1_FEOH', 'HEM-2_EHEM']})
            self.bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
            self.bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

        # temporary mix -> flip charge only for CuA and hemeA

         if specie == 'A_paraccocus' and cycle_state == 'Cua_hemea_R_hemea3_Cub_O':

            ### CHARGE PATCHES
            self.charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
            self.charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
            self.charge_patches.append({'CB4' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
            self.charge_patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})

            ### BOND PATCHES
            self.bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
            self.bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
            self.bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
            self.bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

         if specie == 'A_paraccocus' and cycle_state == 'Cua_hemea_O_hemea3_Cub_R':

            ### CHARGE PATCHES
            self.charge_patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
            self.charge_patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
            self.charge_patches.append({'CB1T' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
            self.charge_patches.append({'A3W2' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})

            ### BOND PATCHES
            self.bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
            self.bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
            self.bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
            self.bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
            self.bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})






