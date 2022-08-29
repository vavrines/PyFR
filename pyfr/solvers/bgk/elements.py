# -*- coding: utf-8 -*-

from ctypes.wintypes import PSIZE
from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np

def setup_BGK(cfg, ndims):
    Nx = cfg.getint('solver', 'Nx', 8)
    Ny = cfg.getint('solver', 'Ny', 8)
    Nz = cfg.getint('solver', 'Nz', 1)

    Lx = cfg.getliteral('solver', 'Lx', [-1, 1])
    Ly = cfg.getliteral('solver', 'Ly', [-1, 1])
    Lz = cfg.getliteral('solver', 'Lz', [-1, 1])


    if ndims == 2:
        nvars = Nx*Ny 
        u = np.zeros((nvars, ndims))
        w = np.zeros((nvars))

        [gauss_pts_x, gauss_wts_x] = np.polynomial.legendre.leggauss(Nx)

        r0 = [0.5*(Lx[0] + Lx[1]), 0.5*(Ly[0] + Ly[1])] 
        rmax = 0.5*(Lx[1] - Lx[0])

        ur = 0.5*(gauss_pts_x + 1)*rmax
        wts_r = (gauss_wts_x/2.0)*ur

        ut = np.linspace(-np.pi, np.pi, Ny, endpoint=False)
        wts_t = np.ones_like(ut)/len(ut)

        ux = r0[0] + np.outer(ur, np.cos(ut))
        uy = r0[1] + np.outer(ur, np.sin(ut))
        wts = np.outer(wts_r, wts_t)            

        u[:,0] = np.reshape(ux, (-1))
        u[:,1] = np.reshape(uy, (-1))
        w = np.reshape(wts, (-1))

        PSint = w*(2*np.pi*rmax)
    elif ndims == 3:
        nvars = Nx*Ny*Nz
        u = np.zeros((nvars, ndims))
        w = np.zeros((nvars))

        [gauss_pts_x, gauss_wts_x] = np.polynomial.legendre.leggauss(Nx)

        r0 = [0.5*(Lx[0] + Lx[1]), 0.5*(Ly[0] + Ly[1]), 0.5*(Lz[0] + Lz[1])] 
        rmax = 0.5*(Lx[1] - Lx[0])

        ur = 0.5*(gauss_pts_x + 1)*rmax
        wts_r = (gauss_wts_x/2.0)*ur
        
        ut = np.linspace(-np.pi, np.pi, Ny, endpoint=False)
        wts_t = np.ones_like(ut)/len(ut)

        # Create open up set
        up = np.linspace(0, np.pi, Nz+1, endpoint=False)
        up = 0.5*(up[1:] + up[:-1])
        wts_p = np.ones_like(up)/len(up)

        ux = r0[0] + np.outer(np.outer(ur, np.sin(up))      , np.cos(ut))
        uy = r0[1] + np.outer(np.outer(ur, np.sin(up))      , np.sin(ut))
        uz = r0[2] + np.outer(np.outer(ur, np.cos(up)), np.ones_like(ut)) 
        wts = np.outer(np.outer(wts_r, wts_t), wts_p)

        u[:,0] = np.reshape(ux, (-1))
        u[:,1] = np.reshape(uy, (-1))
        u[:,2] = np.reshape(uz, (-1))
        w = np.reshape(wts, (-1))

        PSint = w*(2*4./3.*np.pi*rmax**2)

    psi = np.zeros((nvars, ndims+2))
    for i in range(nvars):
        psi[i, 0] = 1
        for j in range(ndims):
            psi[i, 1+j] = u[i, j]
        psi[i, -1] = 0.5*np.linalg.norm(u[i,:])**2
    
    return [u, PSint, psi]
    
def iterate_DVM(Uloc, u, ndims, moments, PSint, gamma, niters):
    def compute_discrete_maxwellian(alpha):
        # Compute the macro/micro velocity defect
        dv2 =  (u[...,0] - alpha[2])**2
        dv2 += (u[...,1] - alpha[3])**2
        if ndims == 3:
            dv2 += (u[...,2] - alpha[4])**2

        M = (alpha[0]*np.exp(-alpha[1]*dv2))
        return M
    
    # Change local variables into alpha vector
    if ndims == 2:
        [rholoc, rhouloc, rhovloc, Eloc] = Uloc
        p = (gamma - 1.)*(Eloc - 0.5*(rhouloc**2 + rhovloc**2)/rholoc)
        theta = p/rholoc
        alpha = np.zeros(ndims+2)

        alpha[0] = rholoc/(2*np.pi*theta)**(ndims/2.0) # A
        alpha[1] = 1.0/(2*theta) # B
        alpha[2] = rhouloc/rholoc # C,D,E
        alpha[3] = rhovloc/rholoc # C,D,E
    elif ndims == 3:
        [rholoc, rhouloc, rhovloc, rhowloc, Eloc] = Uloc
        p = (gamma - 1.)*(Eloc - 0.5*(rhouloc**2 + rhovloc**2 + rhowloc**2)/rholoc)
        theta = p/rholoc
        alpha = np.zeros(ndims+2)

        alpha[0] = rholoc/(2*np.pi*theta)**(ndims/2.0) # A
        alpha[1] = 1.0/(2*theta) # B
        alpha[2] = rhouloc/rholoc # C,D,E
        alpha[3] = rhovloc/rholoc # C,D,E
        alpha[4] = rhowloc/rholoc # C,D,E

    Mloc = compute_discrete_maxwellian(alpha)

    # Perform Newton iterations to find optimal Maxwellian
    for _ in range(niters):
        # Derivatives with respect to alpha
        Q = [None]*(ndims+2)
        Q[0] = 1.0/alpha[0]
        if ndims == 2:
            Q[1] = -((u[...,0] - alpha[2])**2 + (u[...,1] - alpha[3])**2)
            Q[2] = 2*alpha[1]*(u[...,0] - alpha[2])
            Q[3] = 2*alpha[1]*(u[...,1] - alpha[3])
        elif ndims == 3:
            Q[1] = -((u[...,0] - alpha[2])**2 + (u[...,1] - alpha[3])**2 + (u[...,2] - alpha[4])**2)
            Q[2] = 2*alpha[1]*(u[...,0] - alpha[2])
            Q[3] = 2*alpha[1]*(u[...,1] - alpha[3])
            Q[4] = 2*alpha[1]*(u[...,2] - alpha[4])
        
        F = [None]*(ndims+2)
        J = np.zeros((ndims+2, ndims+2))
        for ivar in range(ndims+2):
            psiM = moments[:,ivar]*Mloc
            F[ivar] = np.dot(PSint, psiM) - Uloc[ivar]
            for jvar in range(ndims+2):
                J[ivar, jvar] = np.dot(PSint, Q[jvar]*psiM)

        alpha = alpha - np.linalg.inv(J) @ F
        Mloc = compute_discrete_maxwellian(alpha)
    return Mloc

class BGKElements(BaseAdvectionElements):    
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p'],
                 3: ['rho', 'u', 'v', 'w', 'p']}

    convarmap = {2: ['1'],
                 3: ['1']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p'])]
    }


    def __init__(self, basiscls, eles, cfg):
        self.ndims = eles.shape[2]

        [self.u, self.PSint, self.moments] = setup_BGK(cfg, self.ndims)
        self.nvars = len(self.u)

        super().__init__(basiscls, eles, cfg)

    # Initial Maxwellian state
    def pri_to_con(self, pris, cfg):
        rho, U, p = pris[0], pris[1:-1], pris[-1]
        
        # Multiply velocity components by rho
        rhoUs = [rho*c for c in U]

        # Compute the internal/total energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in U)
        M = np.zeros((self.nupts, self.nvars, self.neles))
        
        # (nupts, _, nelems) = np.shape(M)

        niters = 5 # Large iteration count for ICs
        for uidx in range(self.nupts):
            for eidx in range(self.neles):
                # Get local variables
                rholoc = rho if np.isscalar(rho) else rho[uidx, eidx] 
                rhouloc = rhoUs[0] if np.isscalar(rhoUs[0]) else rhoUs[0][uidx, eidx] 
                rhovloc = rhoUs[1] if np.isscalar(rhoUs[1]) else rhoUs[1][uidx, eidx]
                if self.ndims == 3:
                    rhowloc = rhoUs[2] if np.isscalar(rhoUs[2]) else rhoUs[2][uidx, eidx]
                Eloc = E if np.isscalar(E) else E[uidx, eidx]

                if self.ndims == 2:
                    Uloc = [rholoc, rhouloc, rhovloc, Eloc]
                elif self.ndims == 3:
                    Uloc = [rholoc, rhouloc, rhovloc, rhowloc, Eloc]

                # Compute local Maxwellian
                M[uidx, :, eidx] = iterate_DVM(Uloc, self.u, self.ndims, self.moments, self.PSint, gamma, niters)

        return M

    @staticmethod
    def con_to_pri(cons, cfg, PSint, u, ndims):
        f = cons

        rho = np.dot(PSint, f.swapaxes(0,1))

        Vs = []
        for i in range(ndims):
            rhoU = np.dot(PSint, (f.swapaxes(0,2)*u[:,i]).swapaxes(1,2)).T
            Vs.append(rhoU/rho)

        E = np.dot(PSint, 0.5*(f.swapaxes(0,2)*np.linalg.norm(u, axis=1)**2).swapaxes(1,2)).T

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in Vs))

        return [rho] + Vs + [p]

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.tfluxlin')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.negdivconfbgk')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.limiter')

        ub = self.basis.ubasis
        meanweights = ub.invvdm[:,0]/np.sum(ub.invvdm[:,0])
        self.niters = self.cfg.getint('solver', 'niters')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims, 'nupts': self.nupts, 
            'nvars': self.nvars, 'nverts': len(self.basis.linspts), 
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs, 
            'u': self.u, 'moments': self.moments, 'PSint': self.PSint,
            'srcex': self._src_exprs, 'pi': np.pi,
            'tau': self.cfg.getfloat('constants', 'tau'), 'niters': self.niters,
            'wts': meanweights
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                f=s(self.scal_upts[uin], c), F=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                f=s(self._scal_qpts, c), F=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                f=s(self.scal_upts[uin], l), F=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                f=s(self._scal_qpts, l), F=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )

        self.umat = self._be.const_matrix(self.u)
        self.M = self._be.const_matrix(np.reshape(self.PSint, (1, -1)))
        
        plocsrc = self._ploc_in_src_exprs
        plocupts = self.ploc_at('upts') if plocsrc else None
    
        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconfbgk', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, f=self._scal_upts_cpy,
            u=self.umat, M=self.M
        )

        # Positivity-preserving squeeze limiter
        if self.cfg.getbool('solver', 'limiter', True):
            self.kernels['limiter'] = lambda uin: self._be.kernel(
                'limiter', tplargs=tplargs,
                dims=[self.neles], f=self.scal_upts[uin]
            )