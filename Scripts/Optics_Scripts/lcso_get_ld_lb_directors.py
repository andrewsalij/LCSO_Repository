import berremueller.dielectric_tensor
import numpy as np
import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
import berremueller.berreman_mueller as bm
import berremueller.dielectric_tensor as dt
import ruamel.yaml as yaml
import setup
import os
import Scripts.lcso_params as lcso_params
"Scattering matrix calculations for LCSO w/ various thicknesses. Transfer matrices hit an instability "
filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\TDDFT\ALDA_PBEsol"))
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

dielectric_tensor_init = np.copy(dielectric_tensor)
spec_nm = dt.eV_to_nm(spec)
axis_angle = np.array([1,0,0,-np.pi/2])
quat = dt.axis_angle_to_quaternion(axis_angle)

e_inf = lcso_params.lcso_e_inf
for i in np.arange(3):
    dielectric_tensor[i,i,:] = dielectric_tensor[i,i,:]+e_inf
dielectric_tensor = dt.rotate_2D_tensor(quat,dielectric_tensor,rot_type="quat")

assert dielectric_tensor_init[0,2,:].all() == dielectric_tensor[0,1,:].all(), "rot error"
assert dielectric_tensor_init[1,1,:].all() == dielectric_tensor[2,2,:].all(), "rot error"



wl_to_use = 635
wl_nm_array = np.array([wl_to_use]) #Katherine's laser wavelength
wl_str = str(wl_nm_array[0])

dielectric_near_wl = dielectric_tensor[:,:,np.argmin(np.abs(spec_nm-wl_to_use))]

ld_director, lb_director = yps.get_ld_lb_directors(dielectric_near_wl)
scale_factor = 300
ld_director = ld_director*scale_factor
lb_director = lb_director*scale_factor


