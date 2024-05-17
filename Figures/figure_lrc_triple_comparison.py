import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
from matplotlib.pyplot import  Line2D
import berremueller.dielectric_tensor as dt
import numpy as np
import os
import setup
'''Figure for ACD in inorganics paper--comparison polarization via  ALDA and ALDA+LRC '''

filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\TDDFT\ALDA_PBEsol"))
pol_str_list = ["_X","_Y","_Z"]

data_lrc_pt5 = yps.parse_data(filepath,"LRCpt5",pol_str_list,calc_type="lrc")
data_lrc_1 = yps.parse_data(filepath,"LRC1",pol_str_list,calc_type="lrc")
data_lrc_1pt5 = yps.parse_data(filepath,"LRC1pt5",pol_str_list,calc_type="lrc")

spec_lrc = data_lrc_pt5[0][:,0] #in eV

fname = "lcso_lrc_abs_triple_comparison"

micron_inv_conversion = 1/.1973 #from eV/(hbar c) to micron^-1

spec_lrc_um= dt.eV_to_nm(spec_lrc)*1e-3

abs_lrc_pt5 = yps.convert_to_absorption_coefficient(data_lrc_pt5,spec_lrc,conversion_factor=micron_inv_conversion)
abs_lrc_1 = yps.convert_to_absorption_coefficient(data_lrc_1,spec_lrc,conversion_factor=micron_inv_conversion)
abs_lrc_1pt5 = yps.convert_to_absorption_coefficient(data_lrc_1pt5,spec_lrc,conversion_factor=micron_inv_conversion)


abs_lrc_check = yps.convert_to_absorption_coefficient(data_lrc_pt5,2*np.pi/spec_lrc_um,conversion_factor=1)

#checking unit conversions
assert np.allclose(abs_lrc_pt5[0],abs_lrc_check[0],atol = 1e-3,rtol = 1e-3)

fig,ax = plt.subplots(figsize = (8,3),ncols =3,sharey = True)

abs_lrc_x = [abs_lrc_pt5[0],abs_lrc_1[0],abs_lrc_1pt5[0]]
abs_lrc_y = [abs_lrc_pt5[1],abs_lrc_1[1],abs_lrc_1pt5[1]]
abs_lrc_z = [abs_lrc_pt5[2],abs_lrc_1[2],abs_lrc_1pt5[2]]

colors = ["red","blue","green"]
linestyles = ["solid","--","-."]

base_kwargs = {"linestyles":linestyles,"colors":colors,"legend":False,"ylabel":""}
yps.plot_abs_coef_polarization("",abs_lrc_x,spec_lrc,figure=fig,ax=ax[0],**base_kwargs)
yps.plot_abs_coef_polarization("",abs_lrc_y,spec_lrc,figure=fig,ax=ax[1],**base_kwargs)
yps.plot_abs_coef_polarization("",abs_lrc_z,spec_lrc,figure=fig,ax=ax[2],**base_kwargs)

subtitles = ["a)","b)","c)"]

for i in range(np.size(ax)):
    cur_axis = ax[i]
    cur_axis.set_xlim(0,6)
    cur_axis.set_ylim(-10,100)
    cur_axis.set_xlabel("Energy (eV)")
    cur_axis.text(.04,.93,subtitles[i],transform = cur_axis.transAxes)

ax[0].set_ylabel(r"$\alpha$ ($\mu$m$^{-1}$)")
#legends

plt.rcParams.update({'font.size': 10})

lines = []
for i in range(3):
    lines.append(Line2D([0], [0], linewidth=1, linestyle=linestyles[i],color = colors[i]))
labels = ["0.5","1","1.5"]
leg1 = fig.legend(lines,labels,loc='upper left',bbox_to_anchor=(0.27,.95))
leg1.set_title(r"$\alpha_{xc}$")


fig.subplots_adjust(left=0.12)
fig.savefig(fname+".png",dpi = 500)
fig.show()
