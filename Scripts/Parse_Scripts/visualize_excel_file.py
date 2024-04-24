import xarray

import base_scripts.parse_excel
import numpy as np
import os
import setup
import berremueller.mueller as mueller

mueller_subdir ="Data/Raw_Data/UWM_MuellerMatrix/Excel"
save_dir = "Data/Parsed_Data/UWM_MuellerMatrix/Excel"
mueller_file = "RC2XF_Ka1_MMr_bottom.xlsx"
base_name,ext = os.path.splitext(mueller_file)
save_ext = ".nc"
save_filename = base_name+save_ext
full_filepath = os.sep.join((setup.BASE_DIR,mueller_subdir,mueller_file))
full_savepath = os.sep.join((setup.BASE_DIR,save_dir,save_filename))
os.makedirs(os.sep.join((setup.BASE_DIR,save_dir)),exist_ok = True)

mueller_dataarray = base_scripts.parse_excel.mueller_excel_to_dataarray(full_filepath,verbose=True)

wavelength_label = list(mueller_dataarray.coords.keys())[0]
wavelength_array = mueller_dataarray.coords[wavelength_label].values
mueller_array  =mueller_dataarray.values
wvl_bounds = (400,1500)
idx_start = np.argmin(np.abs(wavelength_array-wvl_bounds[0]))
idx_end = np.argmin(np.abs(wavelength_array-wvl_bounds[1]))
mueller_array = mueller_array[:,:,idx_start:idx_end]
wavelength_array = wavelength_array[idx_start:idx_end]
mueller.mueller_matrix_grid_plot(wavelength_array,mueller_array)


mueller_dataset = mueller_dataarray.to_dataset(name="mueller_matrix")
mueller_dataset.to_netcdf(full_savepath)

load_dataset = xarray.open_dataset(full_savepath)
assert mueller_dataset.equals(load_dataset)