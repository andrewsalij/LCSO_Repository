import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as op
import sklearn.metrics
import os
import setup
import base_scripts.parse_excel as parse_excel

exp_array = parse_excel.load_exp_thickness_numpy()
mean_cd, cd_std_err,mean_thickness,thickness_std_err,front_cd,back_cd,front_std_err,back_std_err\
    = parse_excel.clean_exp_thickness_array(exp_array)
#removing outlier data (thickness > 40 um)
cd_to_use = mean_cd
idx_to_use = np.argwhere(mean_thickness<40).flatten()
cd_to_use = cd_to_use[idx_to_use]
cd_to_use = cd_to_use[~np.isnan(cd_to_use)]
mean_thickness = mean_thickness[idx_to_use]
mean_thickness = mean_thickness[~np.isnan(mean_thickness)]

def get_pol_params(d1, d2, b1, b2):
    d_vec = np.array([d1,d2])
    b_vec = np.array([b1,b2])
    p_vec = b_vec+1j*d_vec
    p_val = np.sqrt(np.dot(p_vec,p_vec))
    i_p = np.imag(p_val)
    r_p = np.real(p_val)
    n_p = np.sqrt(i_p**2+r_p**2)
    return i_p, r_p, n_p
def m03_brown(z,d1, d2, b1, b2,singularity = 1e-8):
    #see (30-31) in https://doi.org/10.1038/s41467-023-44523-1
    # or Brown, C. S. & Bak, A. E. General Lorentz transformation and its
    # application to deriving and evaluating the Mueller matrices of
    # polarization optics. Proc. SPIE 3754, 65–74 (1999).
    i_p,r_p, n_p = get_pol_params(d1, d2, b1, b2)
    a_1 = (np.cosh((i_p*z).astype(float)) - np.cos((r_p*z).astype(float))) / (n_p ** 2+singularity) #polarizance macroscopic parameter
    return 2*32980/np.log(10)*(d1*b2-d2*b1)*a_1

def acd_m03_m00(z,d1, d2, b1, b2):
    i_p,r_p, n_p = get_pol_params(d1, d2, b1, b2)
    d_sq = d1**2+d2**2
    numerator = (d1*b2-d2*b1)*(np.cosh((i_p*z).astype(float)) - np.cos((r_p*z).astype(float)))
    denominator = (r_p**2+d_sq)*np.cosh((i_p*z).astype(float))+(i_p**2-d_sq)*np.cos((r_p*z).astype(float))
    return 2*32980/np.log(10)*numerator/denominator


#pol_array_init = np.array([0.013389717511370684, 0.013202875400325859, 0.017048385109910382, 0.04074044653199602],
                        #  dtype = float) #500 nm wavelength polarizance smm params
pol_array_init = np.array([0.005069503533720284, 0.02163506027617761, -0.0505549674670049, -0.007262965152723405],
                          dtype = float) #635 nm wavelength polarizance smm params
thickness_um_array = np.linspace(0, 40, 101)
thickness_eV_array = thickness_um_array*1e-6/(1.973e-7)
mean_thickness_eV = mean_thickness*1e-6/(1.973e-7)

popt_simple, pcov_simple = op.curve_fit(m03_brown, mean_thickness_eV, cd_to_use, p0=pol_array_init, maxfev=10000)
std_simple= np.sqrt(np.diag(pcov_simple))
cd_fit_simple = m03_brown(thickness_eV_array,*popt_simple)
r2_simple = sklearn.metrics.r2_score(acd_m03_m00(mean_thickness_eV,*popt_simple),cd_to_use)

popt_acd, pcov_acd = op.curve_fit(acd_m03_m00, mean_thickness_eV, cd_to_use, p0=pol_array_init, maxfev=10000,bounds=(-1,np.inf))
std_acd = np.sqrt(np.diag(pcov_acd))
cd_fit_acd = acd_m03_m00(thickness_eV_array,*popt_acd)
r2_acd_full = sklearn.metrics.r2_score(acd_m03_m00(mean_thickness_eV,*popt_acd),cd_to_use)

fig,ax = plt.subplots()
ax.scatter(mean_thickness,cd_to_use,color = "black",label = "Experiment")
ax.plot(thickness_um_array,cd_fit_simple,linestyle = "--",color = "black")
ax.plot(thickness_um_array,cd_fit_acd,linestyle = "dotted",color = "black")
ax.set_xlabel("Thickness (um)")
ax.set_ylabel("CD (mdeg)")
fig.savefig("acd_fig.png")
#fig.show()


save_path = os.sep.join((setup.BASE_DIR,"Data","Modeling_Data","Fit_Thickness"))
os.makedirs(save_path,exist_ok = True)
np.save(os.sep.join((save_path,"thickness_fit_simple.npy")),thickness_um_array)
np.save(os.sep.join((save_path,"cd_fit_simple.npy")),cd_fit_simple)
np.save(os.sep.join((save_path,"thickness_fit_acd.npy")),thickness_um_array)
np.save(os.sep.join((save_path,"cd_fit_acd.npy")),cd_fit_acd)


