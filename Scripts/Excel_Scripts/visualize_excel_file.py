import base_scripts.parse_excel
import matplotlib.pyplot as plt
import os
import setup
import berremueller.mueller as mueller

mueller_subdir ="Data/Raw_Data/UWM_MuellerMatrix/Excel"
mueller_file = "RC2XF_Ka1_MMr_top.xlsx"
full_filepath = os.sep.join((setup.BASE_DIR,mueller_subdir,mueller_file))

mueller_dataarray = base_scripts.parse_excel.mueller_excel_to_dataarray(full_filepath,verbose=True)

wavelength_label = list(mueller_dataarray.coords.keys())[0]
wavelength_array = mueller_dataarray.coords[wavelength_label].values
mueller_array  =mueller_dataarray.values
#normed_mueller_array = mueller.norm_mueller_matrix_stack(mueller_array)
mueller.mueller_matrix_grid_plot(wavelength_array,mueller_array)
