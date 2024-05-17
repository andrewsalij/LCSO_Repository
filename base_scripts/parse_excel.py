import copy
import numpy as np
import pandas as pd
import berremueller.python_util as pu
from typing import Union
import xarray as xr
import re
import warnings
import setup
import os
'''Parse scripts for .xlsx files'''

def load_exp_thickness_numpy(path = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\Thickness_Exp_Data\LCSO_ACD_thickness_data.xlsx"))):
    exp_data = pd.read_excel(path)
    exp_array = exp_data.to_numpy()
    return exp_array

def clean_exp_thickness_array(array):
    '''Note: only for data structured precisely in the same way as LCSO_ACD_thickness_data.xlsx---not a general function
    Returns cd in mdeg, thickness in  um'''
    mean_cd = array[::2,1].astype(float)
    mean_cd = mean_cd[~np.isnan(mean_cd)]
    front_cd = array[::2,3].astype(float)
    back_cd = array[1::2,3].astype(float)
    front_std_err = array[::2,4].astype(float)
    back_std_err = array[1::2,4].astype(float)
    cd_std_err = array[::2,2].astype(float)
    cd_std_err = cd_std_err[~np.isnan(cd_std_err)]
    mean_thickness = array[::2,8].astype(float)
    mean_thickness = mean_thickness[~np.isnan(mean_thickness)]
    thickness_std_err = array[::2,9].astype(float)
    thickness_std_err = thickness_std_err[~np.isnan(thickness_std_err)]
    return mean_cd, cd_std_err,mean_thickness,thickness_std_err,front_cd,back_cd,front_std_err,back_std_err

def parse_mueller_column_label(label:str,idx_start: int = 1):
    mueller_label,remainder_label = split_mueller_column_label(label)
    if mueller_label is not None:
        mueller_indices = extract_digits_from_str(mueller_label)
        return_mueller_indices = tuple(map(lambda x:x-idx_start ,mueller_indices))
    else: return_mueller_indices = (None,None)
    if remainder_label is not None:
        return_remainder_value = extract_number_from_str(remainder_label)
    else: return_remainder_value = None
    return return_mueller_indices,return_remainder_value

def split_mueller_column_label(label: str,idx_start: int = 1)->tuple([Union[str,None],Union[str,None]]):
    all_mueller_matches = re.findall("[mM]_?[\d][\d]",label)
    if len(all_mueller_matches) == 0:
        mueller_match = None
    elif len(all_mueller_matches) == 1:
        mueller_match = all_mueller_matches[0]
    else:
        warnings.warn("Multiple Mueller Matches in label:-early return"+label)
        return None, None
    mueller_label  = None
    if mueller_match is not None:
        mueller_label = mueller_match
    if mueller_label is not None:
        remainder_label = label.replace(mueller_label,"")
    else: remainder_label = None
    if remainder_label == "":remainder_label = None
    return mueller_label,remainder_label

def extract_digits_from_str(string: str)->tuple([int]):
    '''Gets all digits and returns them as a tuple of ints'''
    digits = re.findall('\d',string)
    return tuple(map(int,digits))

def extract_number_from_str(string: str)->int:
    '''Gets first number from string and returns it as a float'''
    number_str = re.match('-?[\d]*(\.[\d]*)?', string).group()
    number = None
    if number_str != "":
        number = float(number_str)
    return number


def mueller_excel_to_dataarray(filepath,verbose = False,normalized = False):
    '''Converts excel file to xarray datarray'''
    dataframe = pd.read_excel(filepath)
    array = dataframe.to_numpy()
    columns = list(dataframe.columns)
    columns_clean = []
    for column in columns:
        clean_column = pu.clean_up_string(column)
        columns_clean.append(clean_column)
    wavelength_idx = pu.get_index_in_list_containing_substring(columns_clean, 'wavelength')
    wavelength_array = copy.deepcopy(array[:,wavelength_idx])
    wavelength_label = columns[wavelength_idx]
    mueller_array = array
    remainder_array = np.zeros(np.size(mueller_array,axis = 1))
    mueller_array_reshaped = np.zeros((4,4,np.size(mueller_array,axis = 0)))
    for i in range(np.size(mueller_array,axis= 1)):
        mueller_label = columns[i]
        mueller_indices,remainder_array[i] = parse_mueller_column_label(mueller_label)
        if mueller_indices[0] is not None:
            if (verbose):
                print("Mueller Indices: ",mueller_indices)
                print("Remainder Value: ",remainder_array[i])
            mueller_array_reshaped[mueller_indices[0],mueller_indices[1],:] = mueller_array[:,i]
    if (normalized):
        mueller_array_reshaped[0,0,:] = 1
    mueller_dataarray = xr.DataArray(mueller_array_reshaped,coords = {wavelength_label:wavelength_array},dims = ['row','column',wavelength_label])
    return mueller_dataarray
