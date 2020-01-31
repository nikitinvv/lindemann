import tomopy
import dxchange
import numpy as np
import h5py
import sys
import os
import json
import tomorectv3d

flat_field_norm = True
flat_field_drift_corr = True 
remove_rings = True
######################################################################################################################

def preprocess_data(prj, flat, dark, FF_norm=flat_field_norm, remove_rings=remove_rings, FF_drift_corr=flat_field_drift_corr, downsapling=0):

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
        prj = tomopy.downsample(prj, level=downsapling)
        prj = tomopy.downsample(prj, level=downsapling, axis=1)
    return prj

if __name__ == "__main__":

    path = sys.argv[1] # path to data folder 
    fid = sys.argv[2] # dataset id, 1127,1128,1136,1137,1149,1150,1156,1157
    proj_start = int(sys.argv[3]) # projection start
    proj_end = int(sys.argv[4]) # projection end # note: ~1005 projections correspond to 180 degrees
    
    #read json file
    with open('info_files.json') as f:
        jfile = json.load(f)
    # read data
    prj, flat, dark, theta = dxchange.read_aps_32id(
        path+jfile[fid]['file_name'], sino=(jfile[fid]['sino_start'], jfile[fid]['sino_end']),proj=(proj_start,proj_end))
    _, flat, _, _ = dxchange.read_aps_32id(# file with flat field
        path+jfile[fid]['file_name_flat'], sino=(jfile[fid]['sino_start'], jfile[fid]['sino_end']),proj=(proj_start,proj_end))        
    # set angles and dark field        
    theta = np.arange(proj_start,proj_start+prj.shape[0])*(theta[1]-theta[0]).astype('float32')
    dark *= 0
    # preprocess
    prj = preprocess_data(prj, flat, dark, FF_norm=flat_field_norm, remove_rings=remove_rings,
                          FF_drift_corr=flat_field_drift_corr, downsapling=jfile[fid]['bin'])                              
    
    # tv recon
    data = np.array(prj.swapaxes(0,1).astype('float32'),order='C')
    data[data<0] = 0
    theta = np.array(theta.astype('float32'),order='C')
    [ns,ntheta,n] = data.shape    
    for k in range(jfile[fid]['center']-10,jfile[fid]['center']+10,2):#1164
        with tomorectv3d.Solver(n, ntheta, ns, jfile[fid]['nsp'], jfile[fid]['method'], jfile[fid]['ngpus'], k/pow(2,jfile[fid]['bin']), jfile[fid]['tv']) as cl:                              
            cl.settheta(theta)
            # reconstruction with 3d tv
            res = np.zeros([ns,n,n],dtype='float32',order='C')
            cl.itertvR(res, data, jfile[fid]['niter'])
            dxchange.write_tiff_stack(res,path+'/rec'+str(fid)+'/'+str(proj_start)+'_'+str(proj_end)+'/rec/res'+str(k)+'.tiff',overwrite=True)    
            with open(path+'/rec'+str(fid)+'/'+str(proj_start)+'_'+str(proj_end)+'/info.json', 'w') as f:
                json.dump(jfile[fid], f)

        