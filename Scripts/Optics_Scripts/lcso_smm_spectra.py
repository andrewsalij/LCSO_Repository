import numpy as np
import matplotlib.pyplot as plt
import base_scripts.yambo_parse_scripts as yps
import berremueller.berreman_mueller as bm
import berremueller.dielectric_tensor as dt
import setup , os
import Scripts.lcso_params as lcso_params
crystal_name_list = ["lcso","ncso","control"]
folder_list = ["ALDA_PBEsol","Na_ALDA_PBEsol","Li_control_ALDA_PBEsol"]
for i in range(len(crystal_name_list)):
    crystal_name = crystal_name_list[i]
    filepath = os.sep.join((setup.BASE_DIR, "Data\Raw_Data\TDDFT", folder_list[i]))
    save_path = os.sep.join((setup.BASE_DIR, r"Data\Modeling_Data\SMM_Set"))
    "Scattering matrix calculations for LCSO. Transfer matrices hit an instability "
    data_set = "alda"  # 'alda' or 'lrc'
    pol_str_list = ["_X", "_Y", "_Z", "_XY", "_XZ"]
    if data_set == "alda":
        folder = "ALDA"
        calc_type = "alda"
    elif data_set == "lrc":
        folder = "LRC1"
        calc_type = "lrc"
    for j in range(2):
        if j == 0: to_flip = False
        else: to_flip = True

        dielectric_tensor,spec = yps.get_3D_dielectric_info_from_yambo_output(filepath,folder,
                                                                              calc_type=calc_type,pol_str_list=pol_str_list,

                                                                   q_str= "1")
        axis_angle = np.array([1,0,0,-np.pi/2])
        quat = dt.axis_angle_to_quaternion(axis_angle)
        if (to_flip):
            quat_2 = dt.axis_angle_to_quaternion(np.array([1,0,0,np.pi]))
            quat = dt.multiply_quaternion(quat,quat_2)
            flip_prefix = "rev"
        else: flip_prefix = ""
        dielectric_tensor = dt.rotate_2D_tensor(quat,dielectric_tensor,rot_type="quat")

        wl_nm_array = dt.eV_to_nm(spec)
        thickness_nm = 300

        e_inf = lcso_params.lcso_e_inf
        for i in np.arange(3):
            dielectric_tensor[i,i,:] = dielectric_tensor[i,i,:]+e_inf

        refl_matrix,trans_matrix = bm.characterize_solo_sample_ps(dielectric_tensor,wl_nm_array,spec,thickness_nm,method = "SM")

        mueller_mat,t_set,t_l,t_r,r_l,r_r = bm.mueller_matrix_suite_from_refl_trans_amplitude_matrices_basis_ps(refl_matrix,trans_matrix)

        cd_mdeg = np.log10(t_r/t_l)*32980

        title = "XZ plane at "+str(thickness_nm)+" nm thickness"
        plt.plot(spec,cd_mdeg)
        plt.ylabel("CD (mdeg)")
        plt.xlabel("Energy (eV)")
        plt.title(title)
        plt.ylim()
        plt.xlim(0,5)
        plt.tight_layout()
        plt.show()

        np.save(os.sep.join((save_path,flip_prefix+crystal_name+"_smm_cd_mdeg.npy")),cd_mdeg)
        np.save(os.sep.join((save_path,flip_prefix+crystal_name+"_smm_spec.npy")),spec)

        abs_avg = (-np.log10(t_r)-np.log10(t_l))/2

        cm_travelled = thickness_nm*1e-7
        abs_coef = abs_avg/cm_travelled

# wavenumber_cm_array = 1/wl_nm_array*1e7
#
# plt.plot(wavenumber_cm_array,abs_coef)
# plt.ylabel("Absorption coef")
# plt.xlabel("Wavenumber (cm^-1)")
# plt.xlim(30000,5000)
# plt.title(title)
# plt.ylim()
# plt.tight_layout()
# plt.show()