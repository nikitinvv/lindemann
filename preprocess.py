import tomopy
import dxchange
import numpy as np
import h5py
import sys
import os
import json
sino_start = 256
sino_end = 1024
flat_field_norm = True
flat_field_drift_corr = True  # Correct the intensity drift
remove_rings = True
binning = 1
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

    path = sys.argv[1]
    fid = sys.argv[2]
    proj_start = int(sys.argv[3])
    proj_end = int(sys.argv[4])
    #read json file
    with open('info_files.json') as f:
        jfile = json.load(f)
    # read data
    prj, flat, dark, theta = dxchange.read_aps_32id(
        path+jfile[fid]['file_name'], sino=(jfile[fid]['sino_start'], jfile[fid]['sino_end']),proj=(proj_start,proj_end))
    _, flat, _, _ = dxchange.read_aps_32id(# file with flat field
        path+jfile[fid]['file_name_flat'], sino=(jfile[fid]['sino_start'], jfile[fid]['sino_end']),proj=(proj_start,proj_end))        
    # set angles and dark field        
    theta = np.arange(proj_start,proj_start+prj.shape[0])*(theta[1]-theta[0])
    dark *= 0
    # preprocess
    prj = preprocess_data(prj, flat, dark, FF_norm=flat_field_norm, remove_rings=remove_rings,
                          FF_drift_corr=flat_field_drift_corr, downsapling=jfile[fid]['bin'])                              
    
    
        
