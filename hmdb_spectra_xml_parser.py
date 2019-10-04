import pandas as pd
from pandas import DataFrame as df
import pickle
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

"""
This script is designed to parse the downloaded spectra .txt files and 
the spectra metadata .xml files from hmdb.

Two dictionaries are output each keyed to unique hmdb id's.  For each 
id, there is a list containing the various spectra as dataframes, or 
the various parsed xml documents as dictionaries.

These will be output, where the ame hmdb key can be used to access the 
molecule metadata from the previously generated pickle.

to do:
1. Okay, so I can parse the spectra and spectra metadata.
    -Handles by each of the two classes  
2. Make two dictionaries with the unique hmdb ID's as keys for each.
3. Have a list in each and append spectra or spectra meta at the same hmdb 
ID as read in.
    -Handled by the third class looping over many of the first two.
    -Necessary as multiple spectra and spectra meta for each hmdb due to 
    repeats at different polarity or with difference CE.
4. Out to pickle.

Next script:
1. Can access the pickle of all of hmdb metadata by hmdb ID.  
2. Chemical formula can be pulled and then exact mass and exact 
mass+adducts can be be calculated.
3. Cross reference from these three pickles togther.
4. First test OH and CO2 losses.
5. Then learn....


"""


class Spectra(object):
    # A class for storing a hmdb .txt spectra as a DataFrame

    def __init__(self, file_name, hmdb_id, spectra):
        self.file_name = file_name
        self.hmdb_id = hmdb_id
        self.spectra_df = spectra_df

    def spec_import(self, file_name):
        spectra_df = pd.read_csv(spectra_path, sep='\t', header=None)
        spectra_df.columns = ['mz', 'i']
        spectra_df['i'] = spectra_df['i'] / spectra_df['i'].max()
        return spectra_df

    def get_hmdb_id(self, file_name):
        pass

    def spec_view(self, spectra_df, file_name):
        x = spectra_df.mz
        y = spectra_df.i
        plt.scatter(x, y)
        plt.title(file_name)
        plt.xlabel('m//z')
        plt.ylabel('intensity')
        plt.show()
        #Save to img?

    def spec_polarity(self):
        # Have to get from other dictionary?  Gross.
        pass

    def spec_has_mz(self):
        # Search with mass tolerance
        pass


class SpecMeta(object):
    # A class for storing a hmdb .xml spectra meta data as a dict

    def __init__(self, file_name, hmdb_id, meta_dict):
        self.file_name = file_name
        self.hmdb_id = hmdb_id
        self.meta_dict = meta_dict
        self.polarity = polarity
        self.highres = highres

    def spec_meta_import(self, file_name):
        tree = ET.parse(file_name)
        root = tree.getroot()
        meta_dict = {}
        for elem in root.iter():
            meta_dict[elem.tag] = elem.text

    def get_hmdb_id(self, file_name, meta_dict):
        pass

    def spec_meta_view(self, meta_dict):
        print(meta_dict)

    def spec_meta_polarity(self, meta_dict):
        return meta_dict['ionization-mode']

    def spec_is_high_res(self, meta_dict):
        # Populate with all headers from db

        low_res = ['Triple_Quad']
        high_res = ['']
        for instrument in low_res:
            if instrument is in meta_dict['notes']:
                return False
        for instrument in high_res:
            if instrument is in meta_dict['notes']:
                return True

class Library(object):
    # A class for a collection of hmdb objects

    def __init__(self, input_folder):
        self.input_folder = input_folder

    def import_loop(self, input_folder):
        # For loop
        pass

    def dual_lookup(self):
        # Calls both objects above
        pass

    def dual_export(self):
        # Exports both libaries to pickle!
        pass


# let's go
# Remember to make an intance of each object from the class...