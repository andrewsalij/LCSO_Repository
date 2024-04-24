import berremueller.dielectric_tensor
import numpy as np
import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
import berremueller.berreman_mueller as bm
import berremueller.dielectric_tensor as dt
import ruamel.yaml as yaml
import setup
import os
"Scattering matrix calculations for LCSO w/ various thicknesses. Transfer matrices hit an instability "
filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\ALDA_PBEsol"))
save_path = os.sep.join((setup.BASE_DIR,r"Data\Modeling_Data\SMM_Thickness"))
os.makedirs(save_path,exist_ok = True)
data_set = "alda"  #'alda' or 'lrc'
pol_str_list = ["_X","_Y","_Z","_XY","_XZ","_YZ"]
if data_set =="alda":
    folder = "ALDA"
    calc_type = "alda"
elif data_set == "lrc":
    folder = "LRC1"
    calc_type = "lrc"

dielectric_tensor,spec = yps.get_3D_dielectric_info_from_yambo_output(filepath,folder,
                                                                      calc_type=calc_type,pol_str_list=pol_str_list,
                                                                      q_str= "1")
axis_angle = np.array([1,0,0,-np.pi/2])
quat = dt.axis_angle_to_quaternion(axis_angle)
axis_angle_2 = np.array([0,0,1,0])
quat2 = dt.axis_angle_to_quaternion(axis_angle_2)
quat = dt.multiply_quaternion(quat,quat2)
axis_angle_3 = np.array([1,0,0,0])
quat = dt.multiply_quaternion(quat,dt.axis_angle_to_quaternion(axis_angle_3))
dielectric_tensor = dt.rotate_2D_tensor(quat,dielectric_tensor,rot_type="quat")


e_inf = 200
for i in np.arange(3):
    dielectric_tensor[i,i,:] = dielectric_tensor[i,i,:]+e_inf

wl_to_use = 500
wl_nm_array = np.array([wl_to_use]) #Katherine's laser wavelength
wl_str = str(wl_nm_array[0])
title_str = str(np.round(quat,decimals=2).tolist())+r"$\epsilon_\infty$ : "+str(e_inf)+" at "+wl_str+" nm"



thickness_nm_array = np.linspace(100,40000,501)
cd_mdeg_array = np.zeros(np.size(thickness_nm_array))
abs_avg_array = np.zeros(np.size(thickness_nm_array))
for i in np.arange(np.size(thickness_nm_array)):
    thickness_nm = thickness_nm_array[i]
    refl_matrix,trans_matrix = bm.characterize_solo_sample_ps(dielectric_tensor,wl_nm_array,spec,thickness_nm)

    mueller_mat,t_set,t_l,t_r,r_l,r_r = bm. mueller_matrix_suite_from_refl_trans_amplitude_matrices_basis_ps(refl_matrix,trans_matrix)
    cd_mdeg = np.log10(t_r/t_l)*32980
    cd_mdeg_array[i] = cd_mdeg
    abs_avg = (-np.log10(t_r)-np.log10(t_l))/2
    abs_avg_array[i] = abs_avg

plt.scatter(thickness_nm_array,cd_mdeg_array)
plt.xlabel("Thickness (nm)")
plt.ylabel("CD (mdeg)")
plt.title(title_str)
plt.savefig(save_path+"smm_cd_lcso_thickness_long.png",dpi = 600)
plt.show()

linear_optics = berremueller.dielectric_tensor.linear_optics_from_dielectric_tensor(dielectric_tensor,spec)

wvl_idx = np.argmin(np.abs(berremueller.dielectric_tensor.eV_to_nm(spec)-wl_to_use))

pol_params = [linear_optics.ld[wvl_idx],linear_optics.ldp[wvl_idx],
              linear_optics.lb[wvl_idx],linear_optics.lbp[wvl_idx]]

np.save(os.sep.join((save_path,"thickness_array.npy")),thickness_nm_array)
np.save(os.sep.join((save_path,"cd_array.npy")),cd_mdeg_array)

yaml_dict = {"data_set":data_set,"wl_nm":wl_to_use,
             "eps_inf_lcso": e_inf,"lcso_quaternion_rotation": quat.tolist(),"lcso_thickness_nm_array": thickness_nm_array.tolist()}
with open(os.sep.join((save_path,"parameters.yaml")), 'w') as file:
    yaml_to_save = yaml.YAML(typ='unsafe',pure=True)
    yaml_to_save.dump(yaml_dict, file)

# plt.scatter(thickness_nm_array,abs_avg_array)
# plt.xlabel("Thickness (nm)")
# plt.ylabel("Abs avg")
# plt.title(title_str)
# plt.savefig("smm_cd_lcso_thickness.png")
# plt.show()