"""
Calculates chemical fingerprints
"""


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs


class FPCalculator(object):
    """
    Calculates Circular, MACCS, topological
    and atom pair fingerprints
    """
    def __init__(self):
        pass

    @staticmethod
    def calculate_cfp(molecular_df, col,
                      bond_radius=2, feature_flag=False):
        """
        Calculates circular fingerprints ECFP/FCFP
        :param molecular_df: pandas data frame containing molecules
        :param col: column with SMILES strings
        :param bond_radius: bond radius for fp generation
        :param feature_flag: FCFP or ECFP
        :return:
        """
        cfp = []

        for index, row in molecular_df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row[col])
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, bond_radius,
                                                           useFeatures=feature_flag)
                cfp.append(fp)
            except:
                cfp.append("N/A")
        if feature_flag:
            molecular_df['fcfp_'+str(bond_radius)] = cfp
        else:
            molecular_df['ecfp_'+str(bond_radius)] = cfp
        return molecular_df

    @staticmethod
    def calculate_topological_fp(molecular_df, col):
        """
        Calculates topological fingperints
        :param molecular_df: pandas data frame containing molecules
        :param col: column with molecules present
        :return:
        """
        fps = []
        for index, row in molecular_df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row[col])
                fp = FingerprintMols.FingerprintMol(mol)
                fps.append(fp)
            except:
                fps.append("N/A")
        molecular_df['top_fp'] = fps
        return molecular_df

    @staticmethod
    def calculate_maccs_key(molecular_df, col):
        """
        Calculates MACCS Key 166 bit fingerprints
        :param molecular_df: pandas data frame containing molecules
        :param col: column with molecules present
        :return:
        """
        fps = []
        for index, row in molecular_df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row[col])
                fp = MACCSkeys.GenMACCSKeys(mol)
                fps.append(fp)
            except:
                fps.append("N/A")
        molecular_df['maccs_166'] = fps
        return molecular_df

    @staticmethod
    def calculate_atom_pair_fp(molecular_df, col):
        """
        Calculates atom pair fingerprint
        :param molecular_df: pandas data frame containing molecules
        :param col: column with molecules present
        :return:
        """
        fps = []
        for index, row in molecular_df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row[col])
                fp = Pairs.GetAtomPairFingerprintAsBitVect(mol)
                fps.append(fp)
            except:
                fps.append('N/A')
        molecular_df['atom_pair_fp'] = fps
        return molecular_df

