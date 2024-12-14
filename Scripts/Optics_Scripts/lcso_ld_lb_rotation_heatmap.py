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
import matplotlib.colors as colors
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

#removing extremely long wavelengths
truncation_idx = 2
dielectric_tensor = dielectric_tensor[:,:,truncation_idx:]
spec = spec[truncation_idx:]
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

angular_resolution = 181
eps_xy = dielectric_tensor[:2, :2,:]  # only calculates for first two indices, in lab_ref frame
angular_set = np.linspace(0, np.pi, angular_resolution)
rotation_matrix = np.zeros((2, 2, angular_resolution))
rotation_matrix[0, 0, :] = np.cos(angular_set)
rotation_matrix[0, 1, :] = -np.sin(angular_set)
rotation_matrix[1, 0, :] = np.sin(angular_set)
rotation_matrix[1, 1, :] = np.cos(angular_set)
eps_xy_rot = np.zeros((2, 2,np.size(eps_xy,axis=-1), angular_resolution), dtype=np.cdouble)
for i in range(angular_resolution):
    eps_xy_rot[:, :,:, i] = np.einsum("ij,jkl->ikl", \
                                    rotation_matrix[:, :, i],
                                    np.einsum("ijl,jk->ikl", eps_xy, np.transpose(rotation_matrix[:, :, i])))
lb_matrix = np.real(np.sqrt(eps_xy_rot[0, 0,:, :]) - np.sqrt(eps_xy_rot[1, 1,:, :]))
ld_matrix = np.imag(np.sqrt(eps_xy_rot[0, 0,:, :]) - np.sqrt(eps_xy_rot[1, 1,:, :]))

angular_set_deg = 180/np.pi*angular_set
x_mesh, y_mesh = np.meshgrid(angular_set_deg,spec_nm)


fig,ax = plt.subplots(nrows= 2,sharex= True)

def gridplot(cur_axis,x_mesh,y_mesh,z_mesh,cmap="bwr",vmin=-1,vmax=1):
    '''domestic--for shared params'''
    gridf= cur_axis.pcolormesh(x_mesh,y_mesh,z_mesh,cmap= cmap,vmin = vmin,vmax = vmax)
    cur_axis.set_ylim(300,1000)
    return gridf
grida = gridplot(ax[0],x_mesh,y_mesh,ld_matrix)
ax[0].set_title("LD")
gridb = gridplot(ax[1],x_mesh,y_mesh,lb_matrix)
ax[1].set_title("LB")

ax[1].set_xlabel("Rotation angle (deg.)")

fig.supylabel("Wavelength (nm)")
norm = colors.Normalize(vmin=-1, vmax=1);
cbar = fig.colorbar(gridb,ax = ax)
cbar.set_label("Signal (arb . u.)")
fig.savefig("ld_lb_lcso_heatmap.png",dpi=1000)

target_wvl = 635 #nm
wvl_idx = np.argmin(np.abs(spec_nm-target_wvl))

fig, ax = plt.subplots()
ax.plot(angular_set_deg,ld_matrix[wvl_idx,:],label = "LD")
ax.plot(angular_set_deg,lb_matrix[wvl_idx,:],label = "LB")
ax.set_ylabel("Signal (arb. u.)")
ax.set_xlabel("Rotation angle (deg.)")
fig.legend(bbox_to_anchor=(0,0,.95,.95))
fig.tight_layout()
fig.savefig("ld_lb_rotation.png",dpi=1000)



