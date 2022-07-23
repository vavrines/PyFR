# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np

def setup_BGK(cfg, ndims):
    Nx = cfg.getint('solver', 'Nx', 8)
    Ny = cfg.getint('solver', 'Ny', 8)
    Nz = cfg.getint('solver', 'Nz', 1)

    quasi1d = (Ny == 1 and Ny == 1)

    Lx = cfg.getliteral('solver', 'Lx', [-1, 1])
    if quasi1d:
        Ly = [-0.5, 0.5]
        Lz = [-0.5, 0.5]
    else:
        Ly = cfg.getliteral('solver', 'Ly', [-0.5, 0.5])
        Lz = cfg.getliteral('solver', 'Lz', [-0.5, 0.5])

    basistype = cfg.get('solver', 'basistype', 'gauss-radial')

    if basistype == 'uniform':
        if ndims == 2:
            nvars = Nx*Ny
            u = np.zeros((nvars, ndims))

            ux = np.linspace(Lx[0], Lx[1], Nx+1)
            ux = 0.5*(ux[1:] + ux[:-1])
            uy = np.linspace(Ly[0], Ly[1], Ny+1)
            uy = 0.5*(uy[1:] + uy[:-1])
            
            uxx, uyy = np.meshgrid(ux, uy)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))
        elif ndims == 3:
            nvars = Nx*Ny*Nz
            u = np.zeros((nvars, ndims))

            ux = np.linspace(Lx[0], Lx[1], Nx+1)
            ux = 0.5*(ux[1:] + ux[:-1])
            uy = np.linspace(Ly[0], Ly[1], Ny+1)
            uy = 0.5*(uy[1:] + uy[:-1])
            uz = np.linspace(Lz[0], Lz[1], Nz+1)
            uz = 0.5*(uz[1:] + uz[:-1])
            
            uxx, uyy, uzz = np.meshgrid(ux, uy, uz)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))
            u[:,2] = np.reshape(uzz, (-1))
    elif basistype == 'gauss-uniform':
        if ndims == 2:
            nvars = Nx*Ny
            u = np.zeros((nvars, ndims))
            w = np.zeros((nvars))

            [gauss_pts_x, gauss_wts_x] = np.polynomial.legendre.leggauss(Nx)
            [gauss_pts_y, gauss_wts_y] = np.polynomial.legendre.leggauss(Ny)

            ux = Lx[0] + 0.5*(Lx[1] - Lx[0])*(gauss_pts_x + 1)
            uy = Ly[0] + 0.5*(Ly[1] - Ly[0])*(gauss_pts_y + 1)
            
            uxx, uyy = np.meshgrid(ux, uy)
            wxx, wyy = np.meshgrid(gauss_wts_x, gauss_wts_y)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))
            w = (1./4.)*np.reshape(wxx, (-1))*np.reshape(wyy, (-1))
        elif ndims == 3:
            nvars = Nx*Ny*Nz
            u = np.zeros((nvars, ndims))
            [gauss_pts_x, gauss_wts_x] = np.polynomial.legendre.leggauss(Nx)
            [gauss_pts_y, gauss_wts_y] = np.polynomial.legendre.leggauss(Ny)
            [gauss_pts_z, gauss_wts_z] = np.polynomial.legendre.leggauss(Nz)

            ux = Lx[0] + 0.5*(Lx[1] - Lx[0])*(gauss_pts_x + 1)
            uy = Ly[0] + 0.5*(Ly[1] - Ly[0])*(gauss_pts_y + 1)
            uz = Lz[0] + 0.5*(Lz[1] - Lz[0])*(gauss_pts_z + 1)
            
            uxx, uyy, uzz = np.meshgrid(ux, uy, uz)
            wxx, wyy, wzz = np.meshgrid(gauss_wts_x, gauss_wts_y, gauss_wts_z)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))
            u[:,2] = np.reshape(uzz, (-1))
            w = (1./8.)*np.reshape(wxx, (-1))*np.reshape(wyy, (-1))*np.reshape(wzz, (-1))
    elif basistype == 'gauss-radial':
        if ndims == 2:
            nvars = Nx*Ny 
            u = np.zeros((nvars, ndims))
            w = np.zeros((nvars))

            [gauss_pts_x, gauss_wts_x] = np.polynomial.legendre.leggauss(Nx)

            r0 = [0.5*(Lx[0] + Lx[1]), 0.5*(Ly[0] + Ly[1])] 
            rmax = 0.5*(Lx[1] - Lx[0])

            ur = 0.5*(gauss_pts_x + 1)*rmax
            ut = np.linspace(-np.pi, np.pi, Ny, endpoint=False)

            wts_t = np.ones_like(ut)/len(ut)
            wts_r = gauss_wts_x*(ur/rmax)

            ux = r0[0] + np.outer(ur, np.cos(ut))
            uy = r0[1] + np.outer(ur, np.sin(ut))
            wts = np.outer(wts_r, wts_t)            

            u[:,0] = np.reshape(ux, (-1))
            u[:,1] = np.reshape(uy, (-1))
            w = np.reshape(wts, (-1))
        elif ndims == 3:
            raise NotImplementedError
    
    if basistype == 'uniform':
        if ndims == 2:
            L = (Lx[1] - Lx[0])*(Ly[1] - Ly[0])
        elif ndims == 3:
            L = (Lx[1] - Lx[0])*(Ly[1] - Ly[0])*(Lz[1] - Lz[0])
        PSint = np.ones(nvars)/nvars*L
    elif basistype == 'gauss-uniform':
        if ndims == 2:
            L = (Lx[1] - Lx[0])*(Ly[1] - Ly[0])
        elif ndims == 3:
            L = (Lx[1] - Lx[0])*(Ly[1] - Ly[0])*(Lz[1] - Lz[0])
        PSint = w*L
    elif basistype == 'gauss-radial':
        r = 0.5*(Lx[1] - Lx[0])
        if ndims == 2:
            A = np.pi*r*r
        # if ndims == 3:
        #     L = (Lx[1] - Lx[0])*(Ly[1] - Ly[0])*(Lz[1] - Lz[0])
        PSint = w*A

    psi = np.zeros((nvars, ndims+2))
    for i in range(nvars):
        psi[i, 0] = 1
        for j in range(ndims):
            psi[i, 1+j] = u[i, j]
        psi[i, -1] = 0.5*np.linalg.norm(u[i,:])**2
    
    return [u, PSint, psi, quasi1d]
    

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

        [self.u, self.PSint, self.moments, self.quasi1d] = setup_BGK(cfg, self.ndims)
        self.K = 1 # Definitely 1 for ndims==2
        self.nvars = len(self.u)

        super().__init__(basiscls, eles, cfg)

    # Initial Maxwellian state
    def pri_to_con(self, pris, cfg):
        rho, U, p = pris[0], pris[1:-1], pris[-1]

        # Compute the macro/micro velocity defect
        dv2 = 0
        for i in range(self.ndims):
            dv2 += (self.u[...,i] - U[i])**2

        # Compute the internal energy
        gamma = cfg.getfloat('constants', 'gamma')
        e = p/(gamma-1)

        if self.quasi1d:
            ndims = 1
        else:
            ndims = self.ndims
        if ndims == 1:
            lam = (0.25)*rho/e
        if ndims == 2:
            lam = (0.5)*rho/e
        if ndims == 3:
            raise NotImplementedError()
            lam = (0.5)*rho/e
        # lam = (0.5)*rho/e

        A = rho*(lam/np.pi)**(ndims/2.0)
        if np.isscalar(A):
            M = A*np.exp(-lam*dv2)
        else:
            M = np.repeat(A[...,np.newaxis], self.nvars, axis=-1)
            for i in range(self.nvars):
                M[..., i] = A*np.exp(-lam*dv2[i])
            M = M.swapaxes(1,2)

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

        tau = self.cfg.getfloat('constants', 'tau')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'u': self.u,
            'moments': self.moments,
            'PSint': self.PSint,
            'K': self.K,
            'srcex': self._src_exprs,
            'pi': np.pi,
            'tau': tau,
            'quasi1d': self.quasi1d
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

        
        plocsrc = self._ploc_in_src_exprs
        plocupts = self.ploc_at('upts') if plocsrc else None
        # solnupts = self._scal_upts_cpy
    
        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconfbgk', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, f=self._scal_upts_cpy
        )