'''
RDK entry point
'''

import pandas as pd
from csds2.chemoinformatics.descriptors.fingerprints import FPCalculator


if __name__ == '__main__':
    df = pd.read_table('/Users/edwardcannon/PycharmProjects/chemogenomics/csds2/data/stahl-cox2.smi',
                       sep=" ", header=0,
                       names=['SMILES'])
    res_df = FPCalculator.calculate_cfp(df,
                                        2,
                                        feature_flag=True)
    print(res_df.head())
    print('topological fps calculation')
    res_df = FPCalculator.calculate_topological_fp(res_df,'SMILES')
    print('now with topological fps')
    print(res_df.head())
    res_df = FPCalculator.calculate_maccs_key(res_df,'SMILES')
    print('maccs')
    print(res_df.head())
    res_df = FPCalculator.calculate_atom_pair_fp(res_df, 'SMILES')
    print(res_df.head())