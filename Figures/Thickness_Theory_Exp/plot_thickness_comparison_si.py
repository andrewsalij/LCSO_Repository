import matplotlib.pyplot as plt
import numpy as np
import xarray
import pandas as pd
import setup
import os
import base_scripts.parse_excel as parse_excel

exp_array = parse_excel.load_exp_thickness_numpy()
mean_cd, cd_std_err,mean_thickness,thickness_std_err,front_cd,back_cd,front_std_err,back_std_err\
    = parse_excel.clean_exp_thickness_array(exp_array)

theory_smm_path = os.sep.join((setup.BASE_DIR,"Data","Modeling_Data","SMM_Thickness"))

theory_fit_path = os.sep.join((setup.BASE_DIR,"Data","Modeling_Data","Fit_Thickness"))

thickness_fit = np.load(os.sep.join((theory_fit_path,"thickness_fit_acd.npy")))
cd_fit = np.load(os.sep.join((theory_fit_path,"cd_fit_acd.npy")))

thickness_theory = np.load(os.sep.join((theory_smm_path,"thickness_array.npy")))
cd_theory = np.load(os.sep.join((theory_smm_path,"cd_array.npy")))


fig,ax = plt.subplots(nrows = 1)

ax.plot(thickness_theory/1e3,cd_theory,color = "blue",label = "SMM")
ax.set_xlabel("Thickness (microns)")
ax.set_ylabel("CD (mdeg)")
ax.set_xlim(0,40)
ax.legend()




fig.subplots_adjust(hspace = 0)
fig.savefig("Figure_3_thickness_c_si")
fig.show()