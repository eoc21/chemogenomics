"""
Calculates chemical descriptors for molecules
"""
from rdkit.Chem import Descriptors
from rdkit import Chem


class MolecularDescriptors(object):
    """
    Base class for molecular descriptor calculation
    """
    def __init__(self):
        pass

    @staticmethod
    def calculate_molecular_descriptors(molecular_df, col):
        """
        Calculates series of molecular descriptors for pandas data frame
        :param molecular_df: pandas data frame
        :param col: molecule column
        :return:
        """
        molecular_df['molecules'] = molecular_df[col].apply(lambda x:
                                                            Chem.MolFromSmiles(x))
        method_list = [func for func in dir(Descriptors)
                       if callable(getattr(Descriptors, func))
                       and not func.startswith("_") and not func.__contains__("version")
                       and not func.__contains__("Fingerprint")]
        for method in method_list:
            molecular_df = MolecularDescriptors.property_calculator(molecular_df,
                                                                    getattr(globals()['Descriptors'], method),
                                                                    str(method))
        return molecular_df

    @staticmethod
    def property_calculator(molecular_df, func_value, col_name):
        """
        Calculates properties of molecules in pandas data frame
        :param molecular_df: pandas data frame
        :param func_value: function to apply
        :param col_name: column name to create
        :return:
        """
        property_values = []
        for index, row in molecular_df.iterrows():
            try:
                property_values.append(func_value(
                    row['molecules']))  # property_values.append(func_value(row['molecules']))
            except:
                property_values.append(0.0)
        molecular_df[col_name] = property_values
        return molecular_df
