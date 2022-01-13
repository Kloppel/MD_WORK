# coding=utf-8
import re

def compare_pkas(f1, f2):

    file1 = open(f1,'r')

    for line in file1:
        if line:
            line = line.strip()
            components = line.split(':')
            resname, resid, segname = re.split('[-_]', components[0])
            pka1 = float(components[1])
            if resname not in ['CTE', 'NTE']:
                identifier = resid + '_' + segname
                file2 = open(f2, 'r')
                for i, line2 in enumerate(file2):
                    if i > 2:
                        line2 = line2.strip()
                        components2 = line2.split(':')
                        # print components2
                        # print len(components2)
                        resname2, resid2, segname2 = re.split('[-_]', components2[0])
                        if resname2 not in ['CTE', 'NTE']:
                            identifier2 = resid2 + '_' + segname2
                            if identifier == identifier2:
                                if '>20' in components2[1]:
                                    components2[1] = 20.0
                                elif '<-10' in components2[1]:
                                    components2[1] = -10.0
                                pka2 = float(components2[1])
                                if abs(pka1-pka2) > 1.5:
                                    print resname + '-' + identifier
                                    # print pka1, pka2
                                    print abs(pka1-pka2)
    return

if __name__ == '__main__':

    file1 = '/scratch/scratch/jdragelj/projects/cco_alexiev/KB_crystal/oxi_nomemb/results.dat'
    # file1 = '/scratch/scratch/jdragelj/projects/cco_alexiev/KB_crystal/oxi_nomemb_tapbs_diff/results.dat'
    file2 = '/media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/awoelke/KB_CcO/cco_paracoccus/3hb3_oxidized_cb4/Results.txt'
    compare_pkas(file1, file2)

