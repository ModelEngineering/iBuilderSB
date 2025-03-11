from iBuilderSB.rhea_maker import RheaMaker # type: ignore

import pandas as pd # type: ignore
import numpy as np
import unittest


IGNORE_TEST = True
IS_PLOT = False


#############################
# Tests
#############################
class TestRheaMaker(unittest.TestCase):

    def setUp(self):
        self.rhea_maker = RheaMaker()

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.rhea_maker.rhea_df, pd.DataFrame))
        self.assertTrue(isinstance(self.rhea_maker.chebi_name_df, pd.DataFrame))

    def testPprintReaction(self):
        if IGNORE_TEST:
            return
        # 10000: Pentanamide + H2O -> pentanoate + NH4+
        result = self.rhea_maker.pprintReaction(10000)
        for element in ["10000", "Pentanamide", "H2O", "pentanoate", "NH4+"]:
            self.assertTrue(element in result)

    def test_stoichiometry_matrix(self):
        #if IGNORE_TEST:
        #    return
        matrix = self.rhea_maker.stoichiometry_matrix
        sym_mat = np.matmul(matrix, matrix.T)
        eigenvalues = np.linalg.eigvals(sym_mat)
        import pdb; pdb.set_trace()

        

if __name__ == '__main__':
    unittest.main()