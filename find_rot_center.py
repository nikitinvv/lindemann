import cupy as cp
import dxchange
import numpy as np
import tomocg as tc
import deformcg as dc
import sys
import os
import matplotlib.pyplot as plt
outdir = './'
def myplot(u, psi, flow):
    [ntheta, nz, n] = psi.shape

    plt.figure(figsize=(20, 14))
    plt.subplot(3, 4, 1)
    plt.imshow(psi[ntheta//4].real, cmap='gray')

    plt.subplot(3, 4, 2)
    plt.imshow(psi[ntheta//2].real, cmap='gray')
    plt.subplot(3, 4, 3)
    plt.imshow(psi[3*ntheta//4].real, cmap='gray')

    plt.subplot(3, 4, 4)
    plt.imshow(psi[-1].real, cmap='gray')

    # plt.subplot(3, 4, 5)
    # plt.imshow(dc.flowvis.flow_to_color(flow[ntheta//4]), cmap='gray')

    # plt.subplot(3, 4, 6)
    # plt.imshow(dc.flowvis.flow_to_color(flow[ntheta//2]), cmap='gray')

    # plt.subplot(3, 4, 7)
    # plt.imshow(dc.flowvis.flow_to_color(flow[3*ntheta//4]), cmap='gray')
    # plt.subplot(3, 4, 8)
    # plt.imshow(dc.flowvis.flow_to_color(flow[-1]), cmap='gray')

    plt.subplot(3, 4, 9)
    plt.imshow(u[nz//2].real)
    plt.subplot(3, 4, 10)
    plt.imshow(u[nz//2+nz//8].real)

    plt.subplot(3, 4, 11)
    plt.imshow(u[:, n//2].real)

    plt.subplot(3, 4, 12)
    plt.imshow(u[:, :, n//2].real)
    if not os.path.exists('tmp'+'_'+str(ntheta)+'/'):
        os.makedirs('tmp'+'_'+str(ntheta)+'/')
    plt.savefig('tmp'+'_'+str(ntheta)+'/flow'+str(k))
    plt.close()
    print(np.linalg.norm(flow))


def update_penalty(psi, h, h0, rho):
    """Update penalty Lagrangian factor rho for faster convergence"""
    r = np.linalg.norm(psi - h)**2
    s = np.linalg.norm(rho*(h-h0))**2
    if (r > 10*s):
        rho *= 2
    elif (s > 10*r):
        rho *= 0.5
    return rho

if __name__ == "__main__":
    name = sys.argv[1]
    part = sys.argv[2]
    #center = np.float(sys.argv[3])
    data = np.load('prj'+name+'_'+part+'.npy').astype('float32')
    theta = np.load('theta'+name+'_'+part+'.npy').astype('float32')
    # Model parameters
    [ntheta,nz,n] = data.shape
    data = data[:,nz//2:nz//2+1]    
    
    # initial guess
    u = np.zeros([1, n, n], dtype='complex64')
    for center in range(n//2-10,n//2+10):
        print(center)
        with tc.SolverTomo(theta, ntheta, 1, n, 1, center) as tslv:
            ucg = tslv.cg_tomo_batch(data, u, 64)
            dxchange.write_tiff_stack(
                ucg.real,  'cg'+'_'+name+'_'+part+'/r'+str(center), overwrite=True)
    