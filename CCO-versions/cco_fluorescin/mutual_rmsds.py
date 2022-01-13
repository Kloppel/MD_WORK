# coding=utf-8

import numpy as np
import cPickle as pickle

class Rmsd(object):

    """This Class is compatible for extracing RMSDS from VMD rmsd output.

       *** all data is in STRING format ***

       Data structure is: {'name': {'structure1': {'structure1' : '0'
                                                   'structure2': 'x'}
                                                                     }
                                                                      }"""

    def __init__(self):

        self.rmsds = {}

    def extract_rmsds(self, name, filename, selection = None):

        #todo: Check for .00 -> make a convention and use it further on!

        if self.check(name):
            raise  AssertionError("This data has been already extracted!")

        self.rmsds[name] = {}
        rmsds_file = open(filename, 'r')


        for line in rmsds_file:
            line = line.split()

            if line == ['ref_mol', 'ref_time', 'mol', 'time', 'rmsd']:
                continue
            if not line:
                continue

            mol_r = line[0]
            mol_c = line[2]

            if mol_r != mol_c:
                print line
                error = "Comparing diferent mol_id-s!"
                raise AssertionError(error)



            frame_r = line[1]
            frame_c = line[3]
            rmsd = line[4]

            # Sometimes VMD produces data with .00 at the end, which is removed in this case
            if '.00' in frame_r:
                frame_r = frame_r[:-3]
            if '.00' in frame_c:
                frame_c = frame_c[:-3]

            if selection is not None:
                if (frame_r not in selection) or (frame_c not in selection):
                    continue

            single_rmsd = {}
            single_rmsd[frame_c] = rmsd

            if frame_r not in self.rmsds[name].keys():
                self.rmsds[name][frame_r] = {}

            self.rmsds[name][frame_r].update(single_rmsd)


    def check(self, name):

        if name in self.rmsds.keys():
            return True
        else:
            return False

    def pickle_results(self, folder, name):

        pickle_results = folder + '/' + name + '.pkl'

        f = open(pickle_results, 'wb')
        pickle.dump(self.rmsds, f)
        f.close()

        return


if __name__ == '__main__':


    filename = '/user/jdragelj/projects/cco_fluorescin/Clustering/data/pyhr_all.dat'

    filter_numerical = np.arange(20,101)
    filter = []
    for frame in filter_numerical:
        filter.append(str(frame))


    test = Rmsd()
    test.extract_rmsds('pyhr', filename, selection = filter)

    print test.rmsds['pyhr']

