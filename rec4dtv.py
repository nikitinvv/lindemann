#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rectv_gpu
import numpy as np
import dxchange
import sys


def getp(a):
    return a.__array_interface__['data'][0]

def takephi(ntheta):
    m = 6  # number of basis functions
    [x, y] = np.meshgrid(np.arange(-ntheta//2, ntheta//2), np.arange(-m//2, m//2))
    phi = np.zeros([m, 2*ntheta], dtype='float32')
    phi[:, ::2] = np.cos(2*np.pi*x*y/ntheta)/np.sqrt(ntheta)
    phi[:, 1::2] = np.sin(2*np.pi*x*y/ntheta)/np.sqrt(ntheta)
    phi[0] = 0  # symmetric
    return phi

if __name__ == "__main__":

    dir = '/data/staff/tomograms/viknik/Lindemann/'
    data = np.load(dir+"bin1prj1149_0.npy").astype('float32')
    theta = np.array(np.load(dir+"bin1theta1149_0.npy").astype('float32'),order='C')

    data = np.array(data.swapaxes(0,1),order='C')
  
    [ns,ntheta,n]=data.shape
    print(data.shape)
    #exit()
    rot_center = 298 # rotation center
    lambda0 = np.float(sys.argv[1])  # regularization parameter 1
    lambda1 = 4  # regularization parameter 2
    nsp = 1 # number of slices to process simultaniously by gpus
    ngpus = 4 # number of gpus
    niter = 256  # number of ADMM iterations
    titer = 4  # number of inner tomography iterations
    
    # take basis functions for decomosition 
    phi = takephi(ntheta) 
    m = phi.shape[0] # number of basis functions
    # creaate class for processing
    cl = rectv_gpu.rectv(n, ntheta, m, ns,
                         nsp, ngpus, rot_center, lambda0, lambda1)
    # angles
    #theta = np.linspace(0, 8*np.pi, ntheta, endpoint=False).astype('float32')  
    # memory for result
    rtv = np.zeros([ns,m,n,n], dtype='float32')
    # Run iterations
    dbg = True # show relative convergence
    cl.run(getp(rtv), getp(data), getp(theta), getp(phi),  niter, titer, dbg)
    # Save result
    for k in range(rtv.shape[0]):
         dxchange.write_tiff_stack(rtv[k], 'rec_tv/rec_'+str(k), overwrite=True)
