import cupy as cp
import dxchange
import numpy as np
import tomocg as tc
import deformcg as dc
import sys
import os
import matplotlib.pyplot as plt
outdir = '/data/staff/tomograms/viknik/lindemann'
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

    data = np.load('prj.npy').astype('float32')#[:,64:96]#[0:6025]
    theta = np.load('theta.npy').astype('float32')#[0:2010]
    # Model parameters
    [ntheta,nz,n] = data.shape
    print(data.shape)
    center = 296 # rotation center
    pnz = 64  # number of slice partitions for simultaneous processing in tomography
    niter = 256
    alpha = 5e-6
    # initial guess
    u = np.zeros([nz, n, n], dtype='complex64')
    psi1 = data.copy()
    lamd1 = np.zeros([ntheta, nz, n], dtype='complex64')
    psi2 = np.zeros([3, nz, n, n], dtype='complex64')    
    lamd2 = np.zeros([3, nz, n, n], dtype='complex64')    
    flow = np.zeros([ntheta, 2], dtype='float32')
    # for center in range(590,610):
    #     print(center)
    #     with tc.SolverTomo(theta, ntheta, nz, n, pnz, center) as tslv:
    #         ucg = tslv.cg_tomo_batch(data, u, 32)
    #         dxchange.write_tiff_stack(
    #             ucg.real,  'cg'+'_'+str(ntheta)+'/rect'+'/r'+str(center), overwrite=True)
    # exit()                
    # ADMM solver
    with tc.SolverTomo(theta, ntheta, nz, n, pnz, center) as tslv:
        # ucg = tslv.cg_tomo_batch(data, u, 64)
        # dxchange.write_tiff_stack(
        #                 ucg.real,  'cg'+'_'+str(ntheta)+'/rect'+'/r', overwrite=True)
        # exit()                        
        with dc.SolverDeform(ntheta, nz, n) as dslv:
            rho1 = 0.5
            rho2 = 0.5
            h01 = psi1
            h02 = psi2
            for k in range(niter):
                # registration
                flow = dslv.registration_shift_batch(data, psi1, 1)     # note swap data and psi for deformation!!!           
                # deformation subproblem
                psi1 = dslv.cg_shift(data, psi1, flow, 2,
                                     tslv.fwd_tomo_batch(u)+lamd1/rho1, rho1)
                #psi1 = (dslv.apply_shift_batch(data, -flow) +
                 #      rho1*tslv.fwd_tomo_batch(u)+lamd1)/(1+rho1)                     
                # regularization subproblem          
                psi2 = tslv.solve_reg(u,lamd2,rho2,alpha)    
                # tomo subproblem
                u = tslv.cg_tomo_batch_ext(psi1-lamd1/rho1, u, 4, rho2/rho1, psi2-lamd2/rho2)

                h1 = tslv.fwd_tomo_batch(u)
                h2 = tslv.fwd_reg(u)
                # lambda update
                lamd1 = lamd1+rho1*(h1-psi1)
                lamd2 = lamd2+rho2*(h2-psi2)
                # checking intermediate results
                #myplot(u, psi, flow)
                if(np.mod(k, 8) == 0):  # check Lagrangian
                    Tpsi = dslv.apply_shift_batch(psi1, flow)
                    lagr = np.zeros(7)
                    lagr[0] = np.linalg.norm(Tpsi-data)**2
                    lagr[1] = np.sum(np.real(np.conj(lamd1)*(h1-psi1)))
                    lagr[2] = rho1*np.linalg.norm(h1-psi1)**2
                    lagr[3] = alpha*np.sum(np.sqrt(np.real(np.sum(psi2*np.conj(psi2), 0))))
                    lagr[4] = np.sum(np.real(np.conj(lamd2*(h2-psi2))))
                    lagr[5] = rho2*np.linalg.norm(h2-psi2)**2
                    lagr[6] = np.sum(lagr[0:5])
                    print(k, np.linalg.norm(flow), rho1, rho2, lagr)
                    # dxchange.write_tiff_stack(
                    #     u.real,  'tmpn'+str(binning)+str(alpha)+'_'+str(ntheta)+'/rect'+str(k)+'/r',overwrite=True)
                    # # dxchange.write_tiff_stack(
                    #     # psi1.real, 'tmp2'+str(binning)+str(alpha)+'_'+str(ntheta)+'/psir'+str(k)+'/r',  overwrite=True)

                    dxchange.write_tiff_stack(
                        u.real,  outdir+'/tmp'+'_'+str(ntheta)+'_'+str(alpha)+'/rect'+str(k)+'/r', overwrite=True)
                    dxchange.write_tiff_stack(
                        psi1.real, outdir+'/tmp'+'_'+str(ntheta)+'_'+str(alpha)+'/psir'+str(k)+'/r',  overwrite=True)
                    np.save(outdir+'/tmp'+'_'+str(ntheta)+'_'+str(alpha)+'/flow'+str(k), flow)                # Updates
                
                rho1 = update_penalty(psi1, h1, h01, rho1)
                rho2 = update_penalty(psi2, h2, h02, rho2)
                h01 = h1.copy()
                h02 = h2.copy()                