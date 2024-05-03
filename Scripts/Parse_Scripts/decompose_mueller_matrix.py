import berremueller.mueller
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

mueller_dataset = xarray.open_dataset(full_savepath)
mueller_dataarray = mueller_dataset["mueller_matrix"]

mueller_decomposed = berremueller.mueller.cloude_decompose_matrix_stack(mueller_dataarray.as_numpy().astype(np.cdouble))

wavelength_label = list(mueller_dataarray.coords.keys())[0]
wavelength_array = mueller_dataarray.coords[wavelength_label].values
mueller_array  =mueller_dataarray.values
wvl_bounds = (400,1500)
idx_start = np.argmin(np.abs(wavelength_array-wvl_bounds[0]))
idx_end = np.argmin(np.abs(wavelength_array-wvl_bounds[1]))
mueller_array_to_plot = mueller_array[:,:,idx_start:idx_end]
wavelength_array_to_plot = wavelength_array[idx_start:idx_end]

mueller.mueller_matrix_grid_plot(wavelength_array_to_plot,mueller_array_to_plot)

mueller.mueller_matrix_grid_plot(wavelength_array_to_plot,mueller_decomposed[:,:,0,idx_start:idx_end]) #plotting first Cloude decomposed matrix
mueller.mueller_matrix_grid_plot(wavelength_array_to_plot,mueller_decomposed[:,:,1,idx_start:idx_end]) #plotting second Cloude decomposed matrix

flip_matrix = np.array([[1,1,1,1],[1,1,1,1],[-1,-1,1,1],[-1,-1,1,1]])
flip_matrix_stack = np.dstack([flip_matrix]*np.size(mueller_decomposed,3))

refl_accounted = mueller_decomposed[:,:,0,:]*flip_matrix_stack
#note-this doesn't have a clean mapping for reflection measurements
logaritmic_decomp = mueller.logm_mueller_matrix_stack(refl_accounted[:,:,idx_start:idx_end])
mueller.mueller_matrix_grid_plot(wavelength_array_to_plot,logaritmic_decomp)