'''Provides access to ChEBI data.'''

from iBuilderSB import constants # type: ignore

import collections
import numpy as np
import pandas as pd # type: ignore
from typing import List, Dict, Tuple

# Columns in ChEBI data
C_REACTION = "reaction"
C_CHEBI = "chebi"
C_SIDE = "side"
C_COMPOUND = "compound"
LEN_COMPOUND = 1 + len(C_COMPOUND)
LEN_CHEBI = 1 + len(C_CHEBI)

# matrix: ndarray
# keys: List[int] - value of the key at the row position
# values: List[int] - value of the value at the column position
OccurrenceMatrixResult = collections.namedtuple("OccurrenceMatrixResult", ["matrix", "keys", "values"])


class RheaMaker(object):

    def __init__(self):
        self.rhea_df = pd.read_csv(constants.RHEA_REACTIONS_PATH)
        self.rhea_df[C_REACTION] = [int(v.split("/")[-1]) for v in self.rhea_df[C_REACTION]]
        self.rhea_df[C_CHEBI] = [int(v.split("/")[-1][LEN_CHEBI:]) for v in self.rhea_df[C_CHEBI]]
        self.rhea_df[C_COMPOUND] = [int(v.split("/")[-1][LEN_COMPOUND:]) for v in self.rhea_df[C_COMPOUND]]
        self.rhea_df[C_SIDE] = [(v.split("/")[-1][-1]) for v in self.rhea_df[C_SIDE]]
        self.reactions = list(set(self.rhea_df[C_REACTION]))
        self.chebis = list(set(self.rhea_df[C_CHEBI]))
        self._stoichiometry_matrix = None
        self.chebi_name_df = pd.read_csv(constants.COMPOUND_NAMES_PATH)

    def pprintReaction(self, rhea_id:int)->str:
        """Prints the reaction with the names of reactants

        Args:
            rhea_id (_type_): _description_

        Returns:
            str: _description_
        """
        ####
        def getChebiName(chebi_id:int)->str:
            names = self.chebi_name_df[self.chebi_name_df["CHEBI_ID"] == chebi_id]["NAME"].values
            return names[0]
        ####
        df = self.rhea_df[self.rhea_df[C_REACTION] == rhea_id]
        lefts = df[df[C_SIDE] == "L"][C_CHEBI].values
        rights = df[df[C_SIDE] == "R"][C_CHEBI].values
        #
        id = str(rhea_id)
        left_names = [str(getChebiName(l)) for l in lefts]
        right_names = [str(getChebiName(l)) for l in rights]
        reactant_str = " + ".join(left_names)
        product_str = " + ".join(right_names)
        return f"{id}: {reactant_str} -> {product_str}"

    @property
    def stoichiometry_matrix(self)->np.ndarray:
        """Creates a stoichiometry matrix for all reactions.
        side="L" is the reactant side of the reaction (negative sign)
        side="R" is the product side of the reaction

        Rank: 10511

        Returns:
            np.ndarray
        """
        if self._stoichiometry_matrix is None:
            self._stoichiometry_matrix = np.zeros((len(self.chebis), len(self.reactions)))
            for chebi in self.chebis:
                chebi_idx = self.chebis.index(chebi)
                df = self.rhea_df[self.rhea_df[C_CHEBI] == chebi]
                for _, row in df.iterrows():
                    reaction_idx = self.reactions.index(row[C_REACTION])
                    sign = 1 if row[C_SIDE] == "R" else -1
                    self._stoichiometry_matrix[chebi_idx, reaction_idx] += sign
        return self._stoichiometry_matrix

    def initialize(self):
        #self.compound_names = pd.read_csv(constants.COMPOUND_NAMES_PATH)
        self.chebi_left_dct = self.makeChebiDct(self.rhea_df[self.rhea_df[C_SIDE] == "L"][C_CHEBI])
        self.chebi_right_dct = self.makeChebiDct(self.rhea_df[self.rhea_df[C_SIDE] == "R"][C_CHEBI])
        self.chebi_both_dct = self.makeChebiDct(self.rhea_df[C_CHEBI])
        import pdb; pdb.set_trace()

    @staticmethod
    def makeOccurrenceMatrix(dct)->OccurrenceMatrixResult:
        """Creates a co-occurrence matrix of keys with values.

        Args:
            dct (dict):
                key: column in matrix
                value: common elements

        Returns:
            np.ndarray: matrix that relates elements
            List[str]: list of elements corresponding to the matrix index
        """
        keys = list(set(dct.keys()))
        values = list(set(dct.values()))
        matrix = np.zeros((len(keys), len(values)))
        for key, value in dct.items():
            matrix[keys.index(key), values.index(value)] += 1
        return OccurrenceMatrixResult(matrix=matrix, keys=keys, values=values)

    def makeChebiDct(self, chebi_ids:List[int])->Dict[int, int]:
        """Makes a dictionary of reactions for the chebi name

        Args:
            chebi_ids (List[str]): _description_

        Returns:
            Dict[str, str]:
                key: chebi_id
                value: rhea reaction id
        """
        chebi_dict:Dict[int, int] = {}
        for chebi_id in chebi_ids:
            chebi_dict[chebi_id] = list(self.rhea_df[self.rhea_df[C_CHEBI] == int(chebi_id)][C_REACTION])   # type: ignore
        return chebi_dict
    
    def calculateStoichiometryMatrixEntropy(self)->pd.Series:
        """Calculates the entropy of the stoichiometry matrix structured as:
        Reactions. Values are entropy.
          "rt_1_0"
          "rt_2_0"
          "rt_0_1"
          "rt_1_1"
          "rt_2_1"
          "rt_0_2"
          "rt_1_2"
          "rt_2_2"
          "rt_other"
        Species dataframe
          <species name>,<chebi_id>, <reactant entropy>, <product entropy>  


        Returns:
            pd.Series: _description_
        """