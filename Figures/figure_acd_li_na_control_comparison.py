import base_scripts.yambo_parse_scripts as yps
import berremueller.dielectric_tensor as dt
import numpy as np
import matplotlib.pyplot as plt
import setup
import os
import Scripts.lcso_params as lcso_params
'''LDLB visualization from quantum espresso and yambo calculations '''
pol_str_list = ["_X","_Z","_XZ"] #order matters
folder = "ALDA"
calc_type = "alda"

filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\TDDFT"))
filepath_li = os.sep.join((filepath,"ALDA_PBEsol"))
filepath_control = os.sep.join((filepath,"Li_control_ALDA_PBEsol"))
filepath_na =os.sep.join((filepath,"Na_ALDA_PBEsol"))

smm_datapath = os.sep.join((setup.BASE_DIR,"Data\Modeling_Data\SMM_Set_100nm"))
lcso_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"lcso_smm_cd_mdeg.npy")))
lcso_rev_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"revlcso_smm_cd_mdeg.npy")))
lcso_smm_spec=  np.load(os.sep.join((smm_datapath,"lcso_smm_spec.npy")))

ncso_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"ncso_smm_cd_mdeg.npy")))
ncso_rev_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"revncso_smm_cd_mdeg.npy")))
ncso_smm_spec=  np.load(os.sep.join((smm_datapath,"ncso_smm_spec.npy")))

control_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"control_smm_cd_mdeg.npy")))
control_rev_smm_cd_mdeg = np.load(os.sep.join((smm_datapath,"revcontrol_smm_cd_mdeg.npy")))
control_smm_spec=  np.load(os.sep.join((smm_datapath,"control_smm_spec.npy")))


data_li = yps.parse_data(filepath_li,folder,pol_str_list,calc_type=calc_type)
data_na = yps.parse_data(filepath_na,folder,pol_str_list,calc_type=calc_type)
data_control = yps.parse_data(filepath_control,folder,pol_str_list,calc_type=calc_type)
spec_li = data_li[0][:,0]
spec_na = data_na[0][:,0]
spec_control = data_control[0][:,0]


def get_both_linear_optics(data,spec_ev,e_inf =lcso_params.lcso_e_inf):
    dielectric_tensor = yps.acd_dat_list_to_2D_dielectric_info(data,e_inf = e_inf)
    linear_optics=  dt.linear_optics_from_dielectric_tensor(dielectric_tensor,spec_ev)
    dielectric_tensor_rev = dt.rotate_2D_tensor(np.array([np.pi,0,0]),dielectric_tensor,rot_type="euler")
    linear_optics_rev = dt.linear_optics_from_dielectric_tensor(dielectric_tensor_rev,spec_ev)
    linear_optics_list = [linear_optics,linear_optics_rev]
    return linear_optics_list

li_linear_optics_list = get_both_linear_optics(data_li,spec_li)
na_linear_optics_list = get_both_linear_optics(data_na,spec_na)
control_linear_optics_list = get_both_linear_optics(data_control,spec_control)

length_nm = 100

fig,axes = plt.subplots(figsize = (6,3),sharey=True,ncols = 2)
kwargs_base = {"figure":fig,"legend_fontsize":9,"ylim":[-2000,2000],"xlim":[0,5],"xlabel":"Energy (eV)","ylabel":""}
yps.plot_acd("",li_linear_optics_list,spec_li,length_nm/197.3,title = "LCSO",ax = axes[0],**kwargs_base)
axes[0].plot(lcso_smm_spec,lcso_smm_cd_mdeg,color = "red",linestyle = "--")
axes[0].plot(lcso_smm_spec,lcso_rev_smm_cd_mdeg,color = "blue",linestyle = "--")

yps.plot_acd("",control_linear_optics_list,spec_control,length_nm/197.3,title = "Control",ax = axes[1],**kwargs_base)
axes[1].plot(control_smm_spec,control_smm_cd_mdeg,color = "red",linestyle = "--")
axes[1].plot(control_smm_spec,control_rev_smm_cd_mdeg,color = "blue",linestyle = "--")

line_front = plt.Line2D([0],[0],color = "red",linestyle = "solid")
line_back = plt.Line2D([0],[0],color = "blue",linestyle = "solid")

line1 = plt.Line2D([0],[0],color = "black",linestyle = "solid")
line2 = plt.Line2D([0],[0],color = "black",linestyle = "--")
for i in range(2):
    legend_init = axes[i].legend([line_front,line_back],["Front","Back"],loc = "upper left")
    legend_add = axes[i].legend([line1,line2],["Pert.","SMM"],loc = "lower left")
    axes[i].add_artist(legend_add)
    axes[i].add_artist(legend_init)
axes[0].set_ylabel("CD (mdeg)",labelpad=-5)
for i in range(2):
    axes[i].text(.04,1.02,["A","B"][i],transform = axes[i].transAxes,fontsize = 10,weight = "bold")
plt.subplots_adjust(left= .2)
fig.savefig("lcso_acd_li_control_comparison_100nm.png",dpi = 500)
fig.show()