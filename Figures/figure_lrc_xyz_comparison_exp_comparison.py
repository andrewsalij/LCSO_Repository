import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
from matplotlib.pyplot import  Line2D
import berremueller.dielectric_tensor as dt
import numpy as np
import os
import setup
'''Figure for ACD in inorganics paper--comparison polarization via  ALDA and ALDA+LRC '''

plot_abs = False
filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\TDDFT\ALDA_PBEsol"))
pol_str_list = ["_X","_Y","_Z"]

data_lrc_pt5 = yps.parse_data(filepath,"LRCpt5",pol_str_list,calc_type="lrc")
data_lrc_1 = yps.parse_data(filepath,"LRC1",pol_str_list,calc_type="lrc")
data_lrc_1pt5 = yps.parse_data(filepath,"LRC1pt5",pol_str_list,calc_type="lrc")
data_alda = yps.parse_data(filepath,"ALDA",pol_str_list,calc_type = "alda")

spec_lrc = data_lrc_pt5[0][:,0] #in eV

if (plot_abs):
    fname = "lcso_lrc_abs_xyz_triple_comparison_eV"
else:
    fname = "lcso_lrc_epsilon_xyz_triple_comparison_eV"

# cm_inv_conversion = 1e4/.1973 #from eV/(hbar c) to cm^-1
#
# spec_lrc_cm_inv= 1/dt.eV_to_nm(spec_lrc)*1e7

convert_factor = 1

abs_lrc_pt5 = yps.convert_to_absorption_coefficient(data_lrc_pt5,spec_lrc,conversion_factor=convert_factor)
abs_lrc_1 = yps.convert_to_absorption_coefficient(data_lrc_1,spec_lrc,conversion_factor=convert_factor)
abs_lrc_1pt5 = yps.convert_to_absorption_coefficient(data_lrc_1pt5,spec_lrc,conversion_factor=convert_factor)
abs_alda = yps.convert_to_absorption_coefficient(data_alda,spec_lrc,conversion_factor=convert_factor)


fig,ax = plt.subplots(figsize = (10,2.5),ncols =4,sharey = True)

colors = ["red","blue","green"]
linestyles = ["solid","--","-."]

base_kwargs = {"linestyles":linestyles,"colors":colors,"legend":False,"xlabel":r"$\omega (eV\hbar^{-1})$","ylabel":""}
if (plot_abs):

    yps.plot_abs_coef_polarization("",abs_alda,spec_lrc,figure = fig,ax = ax[0],**base_kwargs)
    yps.plot_abs_coef_polarization("",abs_lrc_pt5,spec_lrc,figure=fig,ax=ax[1],**base_kwargs)
    yps.plot_abs_coef_polarization("",abs_lrc_1,spec_lrc,figure=fig,ax=ax[2],**base_kwargs)
    yps.plot_abs_coef_polarization("",abs_lrc_1pt5,spec_lrc,figure=fig,ax=ax[3],**base_kwargs)

else:
    yps.plot_dielectric_polarization("",data_alda,spec_lrc,figure = fig,ax = ax[0],**base_kwargs)
    yps.plot_dielectric_polarization("",data_lrc_pt5,spec_lrc,figure = fig,ax = ax[1],**base_kwargs)
    yps.plot_dielectric_polarization("",data_lrc_1,spec_lrc,figure = fig,ax = ax[2],**base_kwargs)
    yps.plot_dielectric_polarization("",data_lrc_1pt5,spec_lrc,figure = fig,ax = ax[3],**base_kwargs)

subtitles = ["A","B","C","D"]

for i in range(np.size(ax)):
    cur_axis = ax[i]
    if (plot_abs):
        cur_axis.set_ylim(0,2e1)
    else:
        cur_axis.set_ylim(0,25)
    cur_axis.text(.04,.93,subtitles[i],transform = cur_axis.transAxes,fontsize= 10, weight=  "bold")
if (plot_abs):
    ax[0].set_ylabel(r"$\alpha (\text{eV}\hbar^{-1}c^{-1})$")
else:
    ax[0].set_ylabel(r"$\epsilon_i$")
#legends

plt.rcParams.update({'font.size': 10})

titles = [r"$\alpha_{xc} = 0$",r"$\alpha_{xc} = 0.5$",r"$\alpha_{xc} = 1.0 $",r"$\alpha_{xc} = 1.5$"]

lines = []
for i in range(3):
    lines.append(Line2D([0], [0], linewidth=1, linestyle=linestyles[i],color = colors[i]))
for i in range(4):
    ax[i].set_xlim(0,5)
    ax[i].set_title(titles[i],pad = -4)
labels = [r"$x$",r"$y$",r"$z$"]
if (plot_abs):
    bbox = (.1,.88)
else:
    bbox = (.13,.88)
leg1 = fig.legend(lines,labels,loc='upper left',bbox_to_anchor=bbox)
leg1.set_title(r"Polarization")
fig.subplots_adjust(left=0.1)
fig.savefig(fname+".png",dpi = 500)
fig.show()

