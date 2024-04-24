import numpy as np
import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
import berremueller.dielectric_tensor as dt
import os
import berremueller.mueller as mueller
import setup
"Scattering matrix calculations for LCSO w/ various thicknesses. Transfer matrices hit an instability "

filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\ALDA_PBEsol"))
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
spec_nm = dt.eV_to_nm(spec)
axis_angle = np.array([1,0,0,np.pi/2])
quat = dt.axis_angle_to_quaternion(axis_angle)
dielectric_tensor = dt.rotate_2D_tensor(quat,dielectric_tensor,rot_type="quat")

e_inf = 100
for i in np.arange(3):
    dielectric_tensor[i,i,:] = dielectric_tensor[i,i,:]+e_inf

wl_nm_to_check = np.array([635]) #Katherine's laser wavelength
wl_idx = np.argmin(np.abs(spec_nm-wl_nm_to_check))

wl_nm_used = spec_nm[wl_idx]
print("Using wavelength: ",wl_nm_used)

lin_opt_params_length1 = dt.linear_optics_from_dielectric_tensor(dielectric_tensor,spec)




polarizance_length1 = mueller.polarizance_from_linear_optics(lin_opt_params_length1)

p_m_array, r_p_array, i_p_array, n_p_array = polarizance_length1.decompose_polarizance()

thickness_nm_array = np.linspace(100,4000,101)
thickness_hceV_array = thickness_nm_array/197.3
b0_array = np.zeros(np.size(thickness_hceV_array))
b1_array = np.zeros(np.size(thickness_hceV_array))
b2_array = np.zeros(np.size(thickness_hceV_array))
b3_array = np.zeros(np.size(thickness_hceV_array))
for i in range(np.size(thickness_hceV_array)):
    thickness_hceV = thickness_hceV_array[i]
    b0,b1,b2,b3 = mueller.brown_params_from_polarizance(polarizance_length1,length = thickness_hceV)
    b0_array[i] = b0[wl_idx]
    b1_array[i]  = b1[wl_idx]
    b2_array[i]  = b2[wl_idx]
    b3_array[i]  = b3[wl_idx]

diff_mat = polarizance_length1.diff_matrix()[:,:,wl_idx]
mean_abs_length1 = lin_opt_params_length1.absorbance
mean_abs_wl = mean_abs_length1[wl_idx]

mueller_mat_obj = mueller.MUELLER_MATRIX(diff_mat)
exp_mat_10um = mueller_mat_obj.macroscopic_matrix(10000/197.3)

ldlb = polarizance_length1.ldlb()
ldlb_at_wl = ldlb[wl_idx]
plt.plot(thickness_nm_array,b0_array,label = "b0")
plt.plot(thickness_nm_array,b1_array,label = "b1")
plt.legend()
plt.show()


d_sum = lin_opt_params_length1.ld[wl_idx]**2+lin_opt_params_length1.ldp[wl_idx]**2
ldlb_v2 = ldlb_at_wl*b1_array/(b0_array+b1_array*d_sum)*32980*2/np.log(10)

limit_value = ldlb_at_wl/(r_p_array[wl_idx]**2+d_sum)*32980*2/np.log(10)

i_p_wl = i_p_array[wl_idx]
r_p_wl = r_p_array[wl_idx]
z_array = thickness_hceV_array
ldlb_v3 = limit_value-limit_value*np.cos(r_p_wl*z_array)/np.cosh(i_p_wl*z_array)
#ldlb_v4 = limit_value-limit_value*np.cos(r_p_wl*z_array)*2*np.exp(i_p_wl*z_array)

plt.plot(thickness_nm_array,ldlb_at_wl*b1_array*32980*2/np.log(10),label = "LDLB")
plt.plot(thickness_nm_array,ldlb_v2,label = "LDLB_v2")
plt.plot(thickness_nm_array,ldlb_v3,label = "LDLB_v3")
#plt.plot(thickness_nm_array,ldlb_v4,label = "LDLB_v4")
plt.plot(thickness_nm_array,limit_value*np.ones(np.size(thickness_nm_array)),linestyle = "dotted",color = "black",label = "Limit")
plt.savefig("ldlb_from_polarizance.png",dpi = 600)
plt.show()

