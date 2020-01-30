import tomopy
import dxchange
import numpy as np
import h5py
import sys

##################################### Inputs #########################################################################
#file_name = '/data/staff/tomograms/viknik/Lindemann/Continuous_1B_8400eV_100ms_1_1149.h5'
#file_name_flat = '/data/staff/tomograms/viknik/Lindemann/Continuous_1B_8400eV_100ms_2_1150.h5'
file_name = '/local/data/2020-01/vnikitin//Continuous_1B_8400eV_100ms_1_1149.h5'
file_name_flat = '/local/data/2020-01/vnikitin//Continuous_1B_8400eV_100ms_2_1150.h5'
sino_start = 300
sino_end = 812
flat_field_norm = True
flat_field_drift_corr = True  # Correct the intensity drift
remove_rings = True
binning = 2
######################################################################################################################


def preprocess_data(prj, flat, dark, FF_norm=flat_field_norm, remove_rings=remove_rings, FF_drift_corr=flat_field_drift_corr, downsapling=binning):

    if FF_norm:  # dark-flat field correction
        prj = tomopy.normalize(prj, flat, dark)
    if FF_drift_corr:  # flat field drift correction
        prj = tomopy.normalize_bg(prj, air=100)
    prj[prj <= 0] = 1  # check dark<data
    prj = tomopy.minus_log(prj)  # -logarithm
    if remove_rings:  # remove rings
        prj = tomopy.remove_stripe_fw(
            prj, level=7, wname='sym16', sigma=1, pad=True)
    if downsapling > 0:  # binning
        prj = tomopy.downsample(prj, level=binning)
        prj = tomopy.downsample(prj, level=binning, axis=1)
    return prj


if __name__ == "__main__":
    proj_start = np.int(sys.argv[1])
    proj_end = proj_start+1005
    # read data
    prj, flat, dark, theta = dxchange.read_aps_32id(
        file_name, sino=(sino_start, sino_end),proj=(proj_start,proj_end))
    _, flat, _, _ = dxchange.read_aps_32id(
        file_name_flat, sino=(sino_start, sino_end),proj=(proj_start,proj_end))        
    # set angles and dark field        
    theta = np.arange(proj_start,proj_start+prj.shape[0])*(theta[1]-theta[0])
    print(theta)
    dark *= 0
    # preprocess
    prj = preprocess_data(prj, flat, dark, FF_norm=flat_field_norm, remove_rings=remove_rings,
                          FF_drift_corr=flat_field_drift_corr, downsapling=binning)                              
    #dxchange.write_tiff(prj,'prj')                          
    np.save('prj'+file_name[-7:-3]+'_'+str(proj_start),prj)        
    np.save('theta'+file_name[-7:-3]+'_'+str(proj_start),theta)  
        
