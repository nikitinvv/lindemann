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
    center = np.float(sys.argv[3])
    data = np.load('prj'+name+'_'+part+'.npy').astype('float32')
    theta = np.load('theta'+name+'_'+part+'.npy').astype('float32')

    [ntheta,nz,n] = data.shape
    pnz = 128  # number of slice partitions for simultaneous processing in tomography
    niter = 128
    # Model parameters
    [ntheta,nz,n] = data.shape
    u = np.zeros([nz, n, n], dtype='complex64')
    psi = data.copy()
    lamd = np.zeros([ntheta, nz, n], dtype='complex64')
    flow = np.zeros([ntheta, 2], dtype='float32')
   
    # ADMM solver
    with tc.SolverTomo(theta, ntheta, nz, n, pnz, center) as tslv:
        with dc.SolverDeform(ntheta, nz, n) as dslv:
            rho = 0.5
            h0 = psi
            for k in range(niter):
                # registration
                flow = dslv.registration_shift_batch(psi, data, 100)
                
                # deformation subproblem
                #psi = dslv.cg_shift(data, psi, flow, 4,
                 #                   tslv.fwd_tomo_batch(u)+lamd/rho, rho)
                psi = (dslv.apply_shift_batch(data, -flow) +
                       rho*tslv.fwd_tomo_batch(u)+lamd)/(1+rho)                                        
                # tomo subproblem
                u = tslv.cg_tomo_batch2(psi-lamd/rho, u, 4)
                h = tslv.fwd_tomo_batch(u)
                # lambda update
                lamd = lamd+rho*(h-psi)

                # checking intermediate results
                if(np.mod(k, 4) == 0):  # check Lagrangian
                    Tpsi = dslv.apply_shift_batch(psi, flow)
                    lagr = np.zeros(4)
                    lagr[0] = np.linalg.norm(Tpsi-data)**2
                    lagr[1] = np.sum(np.real(np.conj(lamd)*(h-psi)))
                    lagr[2] = rho*np.linalg.norm(h-psi)**2
                    lagr[3] = np.sum(lagr[0:3])
                    print(k, np.linalg.norm(flow), rho, lagr)
                    dxchange.write_tiff_stack(
                        u.real,  outdir+'/'+name+'_'+part+'n/tmp'+'_'+str(ntheta)+'_'+'/rect'+str(k)+'/r', overwrite=True)
                    dxchange.write_tiff_stack(
                        psi.real, outdir+'/'+name+'_'+part+'n/tmp'+'_'+str(ntheta)+'_'+'/psir'+str(k)+'/r',  overwrite=True)
                    np.save(outdir+'/'+name+'_'+part+'n/tmp'+'_'+str(ntheta)+'_'+'/flow'+str(k), flow)                # Updates
                rho = update_penalty(psi, h, h0, rho)
                h0 = h
