# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np
from scipy import interpolate
import copy
from pyfr.quadrules import get_quadrule

class BaseFluidElements(object):
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p', 's'],
                 3: ['rho', 'u', 'v', 'w', 'p', 's']}

    convarmap = {2: ['rho', 'rhou', 'rhov', 'E', 's'],
                 3: ['rho', 'rhou', 'rhov', 'rhow', 'E', 's']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p']),
            ('entropy', ['s'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p']),
            ('entropy', ['s'])]
    }

    @staticmethod
    def pri_to_con(pris, cfg, auxvars=False):
        if auxvars:
            return [pris[0]]
        rho, p, s = pris[0], pris[-2], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-2]]

        # Compute the energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[1:-2])

        return [rho] + rhovs + [E] + [s]

    @staticmethod
    def con_to_pri(cons, cfg, auxvars=False):
        if auxvars:
            return [cons[0]]
        rho, E, s = cons[0], cons[-2], cons[-1]

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[1:-2]]

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs))

        return [rho] + vs + [p] + [s]


class EulerElements(BaseFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Register our flux kernel
        self._be.pointwise.register('pyfr.solvers.euler.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.euler.kernels.uflux')

        # Template parameters for the flux kernel
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.cfg.items_as('constants', float))

        if 'flux' in self.antialias:
            	raise ValueError('AA not allowed.')
        else:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts
            )


        self.kernels['tdisu'] = lambda: self._be.kernel(
            'uflux', tplargs=tplargs, dims=[self.nupts, self.neles],
            u=self.scal_upts_inb, smats=self.smat_at('upts'),
            f=self._vect_upts
        )
        # Shock capturing
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')


        if shock_capturing == 'riemann-difference':
            # Register the kernels
            self._be.pointwise.register('pyfr.solvers.euler.kernels.riemanndifference')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.divf_LO')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.subcell')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.residual')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.normalizeresidual')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.limitinterp')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.blendintflux')

            self.kernels['copy_soln_at_fpts'] = lambda: self._be.kernel(
                    'copy', self._scal_fpts_cpy, self._scal_fpts
                )

            self.residual = self._be.matrix((self.nupts, self.neles), tags={'align'}, initval=np.zeros((self.nupts, self.neles)))
            self.alpha_fpts = self._be.matrix((self.nfpts, self.neles), tags={'align'}, initval=np.zeros((self.nfpts, self.neles)))

            rdpts = self.getRDpts()
            diffmatRD = self.generateRDDiffMat(rdpts)
            rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')

            # Template arguments
            tplargs = dict(
                nvars=self.nvars, nupts=self.nupts, nfpts=self.nfpts, ndims=self.ndims,
                c=self.cfg.items_as('constants', float), 
                crd=self.cfg.items_as('solver-riemann-difference', float), 
                order=self.basis.order, diffmatRD=diffmatRD, rsolver=rsolver,
                dt=self.cfg.get('solver-time-integrator', 'dt')
            )

            scpts = self.getSCpts()
            scdiffmat = self.generateSubcellDiffMat(scpts)
            tplargs["scdiffmat"] = scdiffmat
            tplargs["sclengths"] = scpts[1:] - scpts[:-1]

            [self.basis.scptsx, self.basis.scptsy] = self.getAllSCpts()

            self.kernels['divf_LO'] = lambda: self._be.kernel(
                'divf_LO', tplargs=tplargs, dims=[self.neles],
                u=self.scal_upts_inb, plocu=self.ploc_at('upts'), divf=self.scal_upts_outb,
                scsmatsx=self.ele_smat_at('scptsx'), scsmatsy=self.ele_smat_at('scptsy')          
            )


            # self.basis.scpts = self.getSCVertices()
            # self.kernels['divf_LO'] = lambda: self._be.kernel(
            #     'subcell', tplargs=tplargs, dims=[self.neles],
            #     u=self.scal_upts_inb, plocu=self.ploc_at('upts'), divf=self.scal_upts_outb,
            #     scverts=self.ploc_at('scpts'), rcpdjac=self.rcpdjac_at('upts')  
            # )


            self.kernels['residual'] = lambda: self._be.kernel(
                'residual', tplargs=tplargs, dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
                rcpdjac=self.rcpdjac_at('upts'), u=self.scal_upts_inb, r=self.residual
            )

            smoothmat = self.generateSmoothMat()
            tplargs["smoothmat"] = smoothmat
            tplargs["kvars"] = 3
            tplargs["tol"] = 1e-8
            tplargs["quadwts"] = get_quadrule('quad', 'gauss-legendre', (self.basis.order+1)**2).wts

            self.kernels['normalizeresidual'] = lambda: self._be.kernel(
                'normalizeresidual', tplargs=tplargs, dims=[self.neles],
                u=self.scal_upts_inb, r=self.residual, divu=self.scal_upts_outb, 
                usmats=self.ele_smat_at('upts'), rcpdjac=self.rcpdjac_at('upts')
            )

            self.kernels['limitinterp'] = lambda: self._be.kernel(
                'limitinterp', tplargs=tplargs, dims=[self.nfpts, self.neles],
                u=self._scal_fpts
            ) 

            interpmat = self.generateInterpMat(rdpts)
            tplargs["interpmat"] = interpmat
            tplargs["e_max"] = self.cfg.get('solver-riemann-difference', 'e_max')

            [self.basis.rdptsx, self.basis.rdptsy] = self.getAllRDpts()
            tplargs["linesmoothmat"] = self.generate1DSmoothMat()

            tplargs["nrdptsx"] = len(self.basis.rdptsx)
            tplargs["nrdptsy"] = len(self.basis.rdptsx)

            print(self.basis.rdptsx)
            print()
            print(self.basis.rdptsy)

            # Smats order is 
            # [dxi/dx, deta/dx, dxi/dy, deta/dy]
            # [dxi/dx, deta/dx, dzeta/dx, dxi/dy, deta/dy, dzeta/dy,  dxi/dz, deta/dz, dzeta/dz]
            self.kernels['riemanndifference'] = lambda: self._be.kernel(
                'riemanndifference', tplargs=tplargs, dims=[self.neles],
                u=self.scal_upts_inb, plocu=self.ploc_at('upts'), divf=self.scal_upts_outb, res=self.residual, 
                rdsmatsx=self.ele_smat_at('rdptsx'), rdsmatsy=self.ele_smat_at('rdptsy'), alpha_fpts=self.alpha_fpts             
            )

            # Belnds high-order and low-order interface flux based on alpha_fpts
            self.kernels['blendintflux'] = lambda: self._be.kernel(
                'blendintflux', tplargs=tplargs, dims=[self.neles],
                f_HO=self._scal_fpts, f_LO=self._scal_fpts_cpy, alpha_fpts=self.alpha_fpts             
            )

    def getSubcellData(self):
        print(np.shape(self.plocfpts))
        input()

    def getRDpts(self):
        p = self.basis.order
        solpts = self.basis.upts[:p+1,0] # Hack

        rdpts = np.zeros(p+2)
        rdpts[0] = -1.
        rdpts[-1] = 1.
        rdpts[1:-1] = get_quadrule('line', 'gauss-legendre', p).pts

        return rdpts

    def getSCpts(self):
        p = self.basis.order
        solpts = self.basis.upts[:p+1,0] # Hack

        scpts = np.zeros(p+2)
        scpts[0] = -1.
        scpts[-1] = 1.
        scpts[1:-1] = 0.5*(solpts[1:] + solpts[:-1])

        return scpts

    def generateSmoothMat(self):
        p = self.basis.order
        solpts = self.basis.upts[:p+1,0] # Hack
        
        M = np.zeros((self.nupts, self.nupts))

        if self.ndims == 2:
            uidx = lambda xidx, yidx: xidx + yidx*(p+1)
            for i in range(p+1):
                for j in range(p+1):
                    idx = uidx(i,j)
                    ridx = uidx(min(i+1,p),j)
                    lidx = uidx(max(0,i-1),j)
                    tidx = uidx(i,min(j+1,p))
                    bidx = uidx(i,max(0,j-1))

                    dxp = abs(solpts[min(i+1,p)] - solpts[i])
                    dxm = abs(solpts[max(0,i-1)] - solpts[i])
                    dym = abs(solpts[min(j+1,p)] - solpts[j])
                    dyp = abs(solpts[max(0,j-1)] - solpts[j])

                    dxp = 1./dxp if (dxp != 0 and dxm != 0)else 0
                    dxm = 1./dxm if (dxp != 0 and dxm != 0) else 0
                    dyp = 1./dyp if (dyp != 0 and dym != 0)else 0
                    dym = 1./dym if (dyp != 0 and dym != 0) else 0

                    dtot = dxp + dxm + dyp + dym

                    M[idx, ridx] = 0.5*dxp/dtot if dtot != 0 else 0.0
                    M[idx, lidx] = 0.5*dxm/dtot if dtot != 0 else 0.0
                    M[idx, tidx] = 0.5*dyp/dtot if dtot != 0 else 0.0
                    M[idx, bidx] = 0.5*dym/dtot if dtot != 0 else 0.0
                    M[idx, idx] = 0.5 if dtot != 0 else 1.0
        return M

    def generate1DSmoothMat(self):
        p = self.basis.order
        rdpts = self.getRDpts() # Hack
        
        M = np.zeros((p+2, p+2))

        sfac = 0.0
        pfac = 1.0

        if self.ndims == 2:
            # for i in range(1,p+1):
            #     dxp = abs(rdpts[i+1] - rdpts[i])
            #     dxm = abs(rdpts[i-1] - rdpts[i])

            #     dxp = 1./dxp
            #     dxm = 1./dxm 

            #     dtot = dxp + dxm

            #     M[i, i+1] = (1-sfac)*dxp/dtot
            #     M[i, i-1] = (1-sfac)*dxm/dtot
            #     M[i, i] = sfac
            if p < 2:
                M = np.eye(p+2)
            else:
                for i in range(1,p+1):
                    dx = np.abs(rdpts - rdpts[i])**pfac
                    dx[i] = dx[0] = dx[-1] = 1.
                    dx = 1./dx
                    dx[i] = dx[0] = dx[-1] = 0.
                    dtot = np.sum(dx)
                    M[i, :] = (1-sfac)*dx/dtot
                    M[i, i] = sfac
        return M

    def generateInterpMat(self, rdpts):
        p = self.basis.order
        M = np.zeros((p+2, p+1))

        solpts = self.basis.upts[:p+1,0] # Hack
       
        for i in range(p+1):
            vals = np.zeros(p+1)
            vals[i] = 1.
            lag = interpolate.lagrange(solpts, vals)
            for j in range(p+2):
                M[j,i] = lag(rdpts[j])
        return M   

    def generateRDDiffMat(self, rdpts):
        p = self.basis.order
        M = np.zeros((p+1, p+2))

        solpts = self.basis.upts[:p+1,0] # Hack
       
        for i in range(p+2):
            vals = np.zeros(p+2)
            vals[i] = 1.
            dlag = interpolate.lagrange(rdpts, vals).deriv()
            for j in range(p+1):
                M[j,i] = dlag(solpts[j])
        return M

    def generateSubcellDiffMat(self, rdpts):
        p = self.basis.order
        M = np.zeros((p+1, p+2))

        for i in range(p+1):
            dx = rdpts[i+1] - rdpts[i]
            M[i,i] = -1./dx
            M[i,i+1] = 1./dx
        return M

    def generateSolDiffMat(self, rdpts):
        p = self.basis.order
        M = np.zeros((p+1, p+1))

        for i in range(p+1):
            dx = rdpts[i+1] - rdpts[i]
            if i == 0:
                M[i,i] = -1./dx
                M[i,i+1] = 1./dx
            elif i == p:
                M[i,i-1] = -1./dx
                M[i,i] = 1./dx
            else:
                M[i,i-1] = -0.5/dx
                M[i,i+1] = 0.5/dx
        return M

    def getAllRDpts(self):
        p = self.basis.order
        solpts = self.basis.upts[:p+1,0] 
        rdpts = self.getRDpts()

        rdptsx = np.array(np.meshgrid(rdpts, solpts, sparse=False, indexing='xy')).reshape((2,-1)).T
        rdptsy = np.array(np.meshgrid(solpts, rdpts, sparse=False, indexing='xy')).reshape((2,-1)).T

        return [rdptsx, rdptsy]

    def getAllSCpts(self):
        p = self.basis.order
        solpts = self.basis.upts[:p+1,0] 
        rdpts = self.getSCpts()

        rdptsx = np.array(np.meshgrid(rdpts, solpts, sparse=False, indexing='xy')).reshape((2,-1)).T
        rdptsy = np.array(np.meshgrid(solpts, rdpts, sparse=False, indexing='xy')).reshape((2,-1)).T

        return [rdptsx, rdptsy]

    def getSCVertices(self):
        p = self.basis.order
        scpts = self.getSCpts()
        return np.array(np.meshgrid(scpts, scpts, sparse=False, indexing='xy')).reshape((2,-1)).T

