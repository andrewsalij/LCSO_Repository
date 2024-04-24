import numpy as np
import berremueller.berreman_mueller as bm
import berremueller.dielectric_tensor as dt
import base_scripts.yambo_parse_scripts as yps
import berremueller.mueller as mueller
import ruamel.yaml as yaml
import os
import setup

filepath = os.sep.join((setup.BASE_DIR,"Data\Raw_Data\ALDA_PBEsol"))
save_path = os.sep.join((setup.BASE_DIR,r"Data\Modeling_Data\Mueller_Matrix"))
os.makedirs(save_path,exist_ok = True)
data_set = "alda" #'alda' or 'lrc'
pol_str_list = ["_X","_Y","_Z","_XY","_XZ","_YZ"]
if data_set =="alda":
    folder = "ALDA"
    calc_type = "alda"
elif data_set == "lrc":
    folder = "LRC1"
    calc_type = "lrc"

dielectric_tensor,spec_eV = yps.get_3D_dielectric_info_from_yambo_output(filepath,folder,
                                                                      calc_type=calc_type,pol_str_list=pol_str_list,
                                                                      q_str= "1")
axis_angle = np.array([1,0,0,np.pi/2])
quat = dt.axis_angle_to_quaternion(axis_angle)
axis_angle_2 = np.array([0,0,1,.5])
quat2 = dt.axis_angle_to_quaternion(axis_angle_2)
quat = dt.multiply_quaternion(quat,quat2)
dielectric_tensor = dt.rotate_2D_tensor(quat,dielectric_tensor,rot_type="quat")


dielectric_wl_nm_array = dt.eV_to_nm(spec_eV)

e_inf = 150
for i in np.arange(3):
    dielectric_tensor[i,i,:] = dielectric_tensor[i,i,:]+e_inf
force_out_of_plane = False
if (force_out_of_plane):
    dielectric_tensor[2,2,:] = e_inf
    dielectric_tensor[0,2,:] = 0
    dielectric_tensor[1,2,:] = 0
    dielectric_tensor[2,0,:] = 0
    dielectric_tensor[2,1,:] = 0

thickness_nm =  31*1000 #31 microns
def get_jones_dielectric(wl_nm_array,angle_set,dielectric_tensor,dielectric_wl_nm=None,mirror="high_eps"):


    if (dielectric_wl_nm is not None):
        wl_eps_set = dielectric_wl_nm
    else:
        wl_eps_set = wl_nm_array
    lcso_eps = dielectric_tensor

    lcso_thick = thickness_nm
    sample_eps = [lcso_eps]
    sample_thick = [lcso_thick]
    wl_eps_set_list = [wl_eps_set]
    if mirror == "dbr": #puts DBR mirror below sample
        print("Including DBR Mirror:")
        dbr_n_per = 4

        target_lambda_1 = 400
        target_lambda_2 = 600

        e_1, e_2 = 1.5 ** 2, 2.2 ** 2
        dbr_eps_1 = np.identity(3) * e_1
        dbr_eps_2 = np.identity(3) * e_2
        dbr_eps = [dbr_eps_1, dbr_eps_2]

        dbr_wl_nm_set_list = [np.array([0]), np.array([0])]
        dbr_params_obj = bm.DBR_Params_Berre(dbr_n_per, dbr_target_wavelength=target_lambda_1, dbr_eps_list=dbr_eps,
                                             dbr_wl_nm_set_list=dbr_wl_nm_set_list)
        dbr_params_obj_2 = bm.DBR_Params_Berre(dbr_n_per, dbr_target_wavelength=target_lambda_2, dbr_eps_list=dbr_eps,
                                               dbr_wl_nm_set_list=dbr_wl_nm_set_list)

        dbr_params_list = [dbr_params_obj, dbr_params_obj_2]

        sample_params = bm.Sample_Params_Berre(sample_eps, sample_thick, wl_eps_set_list)

        r_matrix_set, t_matrix_set = bm.cavity_berreman_angle_sweep(dbr_params_list, sample_params, wl_nm_array,
                                                                    angle_set, circ=False, n_entry=1, n_exit=1,
                                                                    type="amplitude", half_cavity=True)
        mirror_dict= {"mirror":mirror,"dbr_period":dbr_n_per,"target_lambda_1":target_lambda_1,"target_lambda_2":target_lambda_2,"dbr_eps_1":e_1,"dbr_eps_2":e_2}
    elif mirror == "high_eps": #puts high eps mirror below sample
        print("Including High Epsilon Mirror:")
        high_eps = 10000
        high_eps_tensor = np.identity(3) * high_eps
        mirror_thickess = 10*1000
        sample_thick.extend([mirror_thickess])#10 micron-thick mirror
        sample_eps.extend([high_eps_tensor])
        wl_eps_set_list.extend([np.array([0])])
        r_matrix_set,t_matrix_set = bm.general_berreman_angle_sweep(sample_eps,sample_thick,wl_eps_set_list,wl_nm_array,angle_set,circ=False,matrix_type="ps")
        mirror_dict = {"mirror":mirror,"eps_mirror":high_eps,"mirror_thickness":mirror_thickess}
    else: #no mirror
        r_matrix_set,t_matrix_set = bm.general_berreman_angle_sweep(sample_eps,sample_thick,wl_eps_set_list,wl_nm_array,angle_set,circ=False,matrix_type="ps")
        mirror_dict = {"mirror":mirror}
    return r_matrix_set,t_matrix_set,mirror_dict

wl_nm_array = dielectric_wl_nm_array
angle_set = np.linspace(0, 0, 1) * np.pi / 180  # in radians
angle_set_degrees = angle_set * 180 / np.pi

mirror = "none"
r_matrix_set,t_matrix_set,mirror_dict =  get_jones_dielectric(wl_nm_array,angle_set,dielectric_tensor,dielectric_wl_nm=dielectric_wl_nm_array,mirror=mirror)
if (mirror == "dbr"):
    mirror_prefix = "dbr_"
elif (mirror == "high_eps"):
    mirror_prefix = "high_eps_"
else:
    mirror_prefix = ""
angle_to_check = 0 #degrees
angle_idx = np.argmin(np.abs(angle_set_degrees-angle_to_check))

matrices = [r_matrix_set,t_matrix_set]
matrix_titles = ["reflection","transmission"]
for i in range(2):
    matrix_set = matrices[i]

    jones_matrix_selected = matrix_set[:,:,:,angle_idx]

    mueller_matrix_selected = bm.mueller_from_jones_matrix(jones_matrix_selected)
    mueller_matrix_selected_normed = mueller.norm_mueller_matrix_stack(mueller_matrix_selected)


    file_prefix = mirror_prefix+matrix_titles[i]+"lcso"
    x_bounds = (300,1500)
    mueller.mueller_matrix_grid_plot(wl_nm_array,mueller_matrix_selected,os.sep.join((save_path,file_prefix+str(angle_to_check)+"deg")),x_label = "Wavelength (nm)",mueller_labels="external",x_bounds=x_bounds)
    np.save(os.sep.join((save_path,file_prefix+str(angle_to_check)+"deg")),mueller_matrix_selected_normed)
    mueller.mueller_matrix_grid_plot(wl_nm_array,mueller_matrix_selected_normed,os.sep.join((save_path,file_prefix+"_normed"+str(angle_to_check)+"deg")),
                                     x_label = "Wavelength (nm)",mueller_labels="external",x_bounds=x_bounds,y_bounds = (-1,1))
    np.save(os.sep.join((save_path,file_prefix+"_normed"+str(angle_to_check)+"deg")),mueller_matrix_selected_normed)

yaml_dict = {"eps_inf_lcso": e_inf,"quat_rotation": quat.tolist(),"lcso_thickness_nm": thickness_nm,"x_bounds":x_bounds,"force_out_of_plane_symmetry":force_out_of_plane}
yaml_dict.update(mirror_dict)
with open(os.sep.join((save_path,"parameters.yaml")), 'w') as file:
    yaml_to_write = yaml.YAML(typ='unsafe',pure=True)
    yaml_to_write.dump(yaml_dict, file)