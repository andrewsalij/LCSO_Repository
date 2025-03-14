import yamboparser
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import berremueller.dielectric_tensor as dt

'''Parsing scripts for yambo'''


def parse_yambo_file(filepath,folder,pol_str_list,calc_type="lrc",q_str = "1"):
    f_str_list = []
    for i in range(len(pol_str_list)):
        if calc_type == "alda":
            f_str_list.append("o-ALDA"+pol_str_list[i]+".eps_q"+q_str+"_inv_alda_dyson")
        elif calc_type == "lrc":
            f_str_list.append("o-LRC" + pol_str_list[i] + ".eps_q" + q_str + "_inv_lrc_dyson")
    files = []
    for i in range(len(f_str_list)):
        full_path = os.sep.join((filepath,folder))
        files.append(yamboparser.YamboFile(f_str_list[i],full_path,**{"zip_tags" : True}))
    return files
def load_data(file):
    data = file.data
    data = pd.DataFrame(data).to_numpy()
    return data
def load_data_list(file_list):
    data_list = []
    for i in range(len(file_list)):
        data_list.append(load_data(file_list[i]))
    return data_list
def parse_data(filepath,folder,pol_str_list,calc_type,q_str = "1"):
    files = parse_yambo_file(filepath,folder,pol_str_list,calc_type,q_str)
    data_list = load_data_list(files)
    return data_list

def check_spec_constistency(data_list):
    spec = data_list[0][:,0]
    for data in data_list:
        if np.array_equal(data[:,0],spec) == False:
            raise ValueError("Spectra are not the same")

def acd_dat_list_to_2D_dielectric_info(data_list,e_inf = 0):
    '''x and z indices used as default--can be any two indices (x,y,z)'''
    check_spec_constistency(data_list)
    eps_xx = data_list[0][:,2]+1j*data_list[0][:,1]
    eps_zz = data_list[1][:,2]+1j*data_list[1][:,1]
    eps_xz = (data_list[2][:,2]+1j*data_list[2][:,1]-(eps_xx+eps_zz)/2)
    dielectric_tensor = np.zeros((3,3,np.size(eps_xz)),dtype = np.cdouble)
    dielectric_tensor[0,0,:],dielectric_tensor[1,1,:]= eps_xx+e_inf,eps_zz+e_inf
    dielectric_tensor[0,1,:],dielectric_tensor[1,0] = eps_xz,eps_xz
    dielectric_tensor[2,2,:] = 1
    return dielectric_tensor

def get_3D_eps_indices(pol_str_list):
    xx_index,yy_index,zz_index,xy_index,xz_index,yz_index = None,None,None,None,None,None
    for i in range(len(pol_str_list)):
        pol_str=  pol_str_list[i].lower()
        if ("x" in pol_str):
            if ("y" in pol_str):
                xy_index = i
            elif ("z" in pol_str):
                xz_index = i
            else:
                xx_index = i
        elif ("y" in pol_str):
            if ("z" in pol_str):
                yz_index = i
            else:
                yy_index = i
        elif ("z" in pol_str):
            zz_index = i
    indices_to_return=[xx_index,yy_index,zz_index,xy_index,xz_index,yz_index]
    for i in range(5):
        if i not in indices_to_return:
            raise ValueError("index missing from 3D epsilon list")
    return indices_to_return

def full_data_list_to_3D_dielectric_info(data_list,pol_str_list):
    '''Returns full dielectric tensor and spectrum.
    Input data in order x,y,z, xy, xz, yz'''
    check_spec_constistency(data_list)
    eps_indices = get_3D_eps_indices(pol_str_list) #fixes if indices are out of order
    spec_ev = data_list[0][:,0]
    if (eps_indices[0] is not None):
        eps_xx = data_list[eps_indices[0]][:,2]+1j*data_list[eps_indices[0]][:,1]
    else: eps_xx = np.ones(np.size(spec_ev),dtype = np.cdouble)
    if (eps_indices[1] is not None):
        eps_yy = data_list[eps_indices[1]][:,2]+1j*data_list[eps_indices[1]][:,1]
    else: eps_yy = np.ones(np.size(spec_ev),dtype = np.cdouble)
    if (eps_indices[2] is not None):
        eps_zz = data_list[eps_indices[2]][:,2]+1j*data_list[eps_indices[2]][:,1]
    else: eps_zz = np.ones(np.size(spec_ev),dtype = np.cdouble)
    if (eps_indices[3] is not None):
        eps_xy = (data_list[eps_indices[3]][:,2]+1j*data_list[eps_indices[3]][:,1])-(eps_xx+eps_yy)/2
    else: eps_xy = np.zeros(np.size(spec_ev),dtype = np.cdouble)
    if (eps_indices[4] is not None):
        eps_xz = (data_list[eps_indices[4]][:,2]+1j*data_list[eps_indices[4]][:,1])-(eps_xx+eps_zz)/2
    else: eps_xz = np.zeros(np.size(spec_ev),dtype = np.cdouble)
    if (eps_indices[5] is not None):
        eps_yz = (data_list[eps_indices[5]][:,2]+1j*data_list[eps_indices[5]][:,1])-(eps_yy+eps_zz)/2
    else: eps_yz = np.zeros(np.size(spec_ev),dtype = np.cdouble)
    dielectric_tensor = np.zeros((3,3,np.size(spec_ev)),dtype = np.cdouble)
    dielectric_tensor[0,0,:],dielectric_tensor[1,1,:],dielectric_tensor[2,2:] = eps_xx,eps_yy,eps_zz
    dielectric_tensor[0,1,:],dielectric_tensor[1,0] = eps_xy,eps_xy
    dielectric_tensor[0,2,:],dielectric_tensor[2,0,:] = eps_xz,eps_xz
    dielectric_tensor[1,2,:],dielectric_tensor[2,1,:] = eps_yz,eps_yz
    return dielectric_tensor,spec_ev
def finish_plot(ax,filename,figure,xlim,ylim,xlabel,ylabel,title,legend,legend_fontsize,tight_layout):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if tight_layout == True: figure.tight_layout()
    if legend == True: legend_init = ax.legend(fontsize=legend_fontsize)
    figure.savefig(filename,dpi= 500)
    figure.show()

def plot_acd(filename,linear_optics_list,spec,length,
             labels = ["Front","Back"],colors = ["red","blue"],
             unit_type = "mdeg",
             xlim= (40000,5000),ylim =(-20000,20000),
             xlabel= r"Wavenumber (cm$^{-1}$)",ylabel= "ACD (mdeg)",
             title="",legend=True,legend_fontsize=16,tight_layout=True,
             figure= None,ax=None,detail_level=1,linestyles=  ["solid","solid"]):
    if (unit_type == "mdeg"):
        unit_factor = 32980*2/np.log(10)
    else:
        unit_factor = 1
    #TODO add in other unit types
    if detail_level ==1:
        length_factor = length**2 #acd normally-1/2 factor is in linear_optics()
    elif detail_level ==2:
        #TODO add in factor for higher order Mueller terms
        length_factor = 0
    else:raise ValueError("Invalid detail level")
    total_length_factor = unit_factor*length_factor
    if figure == None:figure = plt.figure()
    if ax == None:ax = figure.add_subplot()
    for i in range(len(linear_optics_list)):
        ax.plot(spec,total_length_factor*linear_optics_list[i].ldlb(),label=labels[i],color = colors[i],linestyle = linestyles[i])
    finish_plot(ax,filename,figure,xlim,ylim,xlabel,ylabel,title,legend,legend_fontsize,tight_layout)

def plot_dielectric_polarization(filename,data_list,spec,labels = None,
                                 colors = ["red","green","blue"],linestyles = ["solid","solid","solid"],
                                 xlim= (50000,1000),ylim =(0,20),
                                 xlabel= r"Wavenumber (cm$^{-1}$)",ylabel= r"$\epsilon_i$",
                                 title="",legend=True,legend_fontsize=16,tight_layout=True,
                                 figure= None,ax=None):
    if figure == None:figure = plt.figure()
    if ax == None:ax = figure.add_subplot()
    if labels == None:labels = ["" for i in range(len(data_list))]
    for i in range(len(data_list)):
        ax.plot(spec,data_list[i][:,1],label=labels[i],color = colors[i],linestyle= linestyles[i])
    finish_plot(ax,filename,figure,xlim,ylim,xlabel,ylabel,title,legend,legend_fontsize,tight_layout)

def plot_abs_coef_polarization(filename,abs_coef_list,spec,labels = None,
                                 colors = ["red","green","blue"],linestyles = ["solid","solid","solid"],
                                 xlim= (50000,1000),ylim =(0,20),
                                 xlabel= r"Wavenumber (cm$^{-1}$)",ylabel= r"$\alpha_i  (cm$^{-1}$)$",
                                 title="",legend=True,legend_fontsize=16,tight_layout=True,
                                 figure= None,ax=None):
    if figure == None:figure = plt.figure()
    if ax == None:ax = figure.add_subplot()
    if labels == None:labels = ["" for i in range(len(abs_coef_list))]
    for i in range(len(abs_coef_list)):
        ax.plot(spec,abs_coef_list[i],label=labels[i],color = colors[i],linestyle= linestyles[i])
    finish_plot(ax,filename,figure,xlim,ylim,xlabel,ylabel,title,legend,legend_fontsize,tight_layout)
def convert_to_absorption_coefficient(data_list,spec,conversion_factor = 1,force_eps_r_1=False):
    abs_coef_list = []
    for i in range(len(data_list)):
        if force_eps_r_1 == True:
            eps = 1+data_list[i][:,1]*1j
        else:
            eps = data_list[i][:, 2] + data_list[i][:, 1]*1j
        ref_index = np.sqrt(eps)
        abs_coef = 2 * np.imag(ref_index) * spec / dt.unit_defs_base.c  # in units of inverse length
        abs_coef = conversion_factor*abs_coef
        abs_coef_list.append(abs_coef)
    return abs_coef_list
def data_prime_reconstructed(data_list,data_rot_list):
    '''
    takes data and data_rot, which are of the form A,B,AB and A',B',A'B' where ' is at a 45 degree counterclockwise rotation
    from the initial and reorders to A,B,A', A',B',A, which should be in principle nearly identical to the original data objects
    :return:
    '''
    data_a,data_b,data_ab = data_list[0],data_list[1],data_list[2]
    data_a_rot,data_b_rot,data_ab_rot = data_rot_list[0],data_rot_list[1],data_rot_list[2]
    data_abap = [data_a,data_b,data_a_rot]
    data_apbpa= [data_a_rot,data_b_rot,data_a]
    return data_abap,data_apbpa

def get_3D_dielectric_info_from_yambo_output(filepath,folder,
                                          pol_str_list = ["_X","_Y","_Z","_XY","_XZ","_YZ"],
                                          calc_type="alda",q_str = "1"):
    data_list = parse_data(filepath,folder,pol_str_list,calc_type,q_str)
    dielectric_info,spec_ev = full_data_list_to_3D_dielectric_info(data_list,pol_str_list)
    return dielectric_info,spec_ev

def angle_to_director(angle):
    return np.array([np.cos(angle),np.sin(angle),0])
def get_ld_lb_directors(dielectric_tensor,angular_resolution = 180):
    '''For a dielectric tensor, returns the vectors that maximize LD and LB for the 0,1 axes'''
    eps_xy = dielectric_tensor[:2,:2] # only calculates for first two indices
    angular_set = np.linspace(0,np.pi,angular_resolution)
    rotation_matrix = np.zeros((2,2,angular_resolution))
    rotation_matrix[0,0,:] = np.cos(angular_set)
    rotation_matrix[0,1,:] = -np.sin(angular_set)
    rotation_matrix[1,0,:] = np.sin(angular_set)
    rotation_matrix[1,1,:] = np.cos(angular_set)
    eps_xy_rot = np.zeros((2,2,angular_resolution),dtype = np.cdouble)
    for i in range(angular_resolution):
        eps_xy_rot[:,:,i] = np.einsum("ij,jk->ik", \
              rotation_matrix[:,:,i], np.einsum("ij,jk->ik", eps_xy, np.transpose(rotation_matrix[:,:,i])))
    lb_array = np.real(np.sqrt(eps_xy_rot[0,0,:])-np.sqrt(eps_xy_rot[1,1,:]))
    ld_array = np.imag(np.sqrt(eps_xy_rot[0,0,:])-np.sqrt(eps_xy_rot[1,1,:]))
    plt.plot(angular_set,ld_array,label = "LD")
    plt.plot(angular_set,lb_array,label = "LB")
    plt.xlabel("Rotation angle (rad)")
    plt.ylabel("LD/LB (a.u.)")
    plt.legend()
    plt.savefig("ld_lb_rotation.png",dpi=500)
    plt.show()
    ld_max_idx = np.argmax(ld_array)
    lb_max_idx = np.argmax(lb_array)
    ld_max_angle = angular_set[ld_max_idx]
    lb_max_angle = angular_set[lb_max_idx]
    ld_max_size = np.abs(ld_array[ld_max_idx])
    lb_max_size = np.abs(lb_array[lb_max_idx])
    return angle_to_director(ld_max_angle)*ld_max_size,angle_to_director(lb_max_angle)*lb_max_size


