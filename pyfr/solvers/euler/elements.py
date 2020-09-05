# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np
from scipy import interpolate
import copy
from pyfr.quadrules import get_quadrule

class BaseFluidElements(object):
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p'],
                 3: ['rho', 'u', 'v', 'w', 'p']}

    convarmap = {2: ['rho', 'rhou', 'rhov', 'E'],
                 3: ['rho', 'rhou', 'rhov', 'rhow', 'E']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p'])]
    }

    @staticmethod
    def pri_to_con(pris, cfg, auxvars=False):
        if auxvars:
            return [pris[0]]
        rho, p = pris[0], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-1]]

        # Compute the energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[1:-1])

        return [rho] + rhovs + [E]

    @staticmethod
    def con_to_pri(cons, cfg, auxvars=False):
        if auxvars:
            return [cons[0]]
        rho, E = cons[0], cons[-1]

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[1:-1]]

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs))

        return [rho] + vs + [p]


class EulerElements(BaseFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Register our flux kernel
        self._be.pointwise.register('pyfr.solvers.euler.kernels.tflux')

        # Template parameters for the flux kernel
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.cfg.items_as('constants', float))

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts
            )
        else:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts
            )

        # Shock capturing
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')

        # Shockcell is an elementwise value (0,1) dictating if there is a shock in the cell (1)
        self.shockcell = self._be.matrix((1, self.neles), tags={'align'}, initval=np.zeros((1, self.neles)))

        if shock_capturing == 'riemann-difference':
            # Register the kernels
            self._be.pointwise.register('pyfr.solvers.euler.kernels.rdshocksensor')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.riemanndifference')

            self.kernels['copy_soln_at_fpts'] = lambda: self._be.kernel(
                    'copy', self._scal_fpts_cpy, self._scal_fpts
                )

            # Sense shocks using density 
            shockvar = 0
            shocksensor = self.cfg.get('solver-riemann-difference', 'shock-sensor', 'none')

            # Get matrices for convergence sensor
            pstages = self.cfg.getliteral('solver-riemann-difference', 'order-stages', [1,2,3])
            if max(pstages) > self.basis.order:
                raise ValueError('Invalid projection order for convergence stages. {0}'.format(pstages))

            [diffxi, diffeta, diffzeta, proj1, proj2, proj3, quadwts] = self.generateProjectionMats(pstages)


            # Sense shocks and choose FR or RD
            self.kernels['rdshocksensor'] = lambda: self._be.kernel(
               'rdshocksensor', tplargs=tplargs, dims=[self.neles], u=self.scal_upts_inb, 
               shockcell=self.shockcell, divf_fr=self.scal_upts_outb, divf_rd=self._scal_upts_cpy,
               usmats=self.ele_smat_at('upts')
            )


            diffmat = self.generateDiffMat()
            rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')

            # Template arguments
            tplargs = dict(
                nvars=self.nvars, nupts=self.nupts, nfpts=self.nfpts, ndims=self.ndims,
                c=self.cfg.items_as('constants', float), shocksensor=shocksensor,
                crd=self.cfg.items_as('solver-riemann-difference', float), 
                order=self.basis.order, svar=shockvar, diffmat=diffmat, rsolver=rsolver,
                vdm=self.basis.ubasis.vdm.T, 
                invvdm=self.basis.ubasis.invvdm.T, 
                ubdegs=[sum(dd) for dd in self.basis.ubasis.degrees],
                l1degs=[max(dd) for dd in self.basis.ubasis.degrees],
                diffxi=diffxi, diffeta=diffeta, diffzeta=diffzeta,
                proj1=proj1, proj2=proj2, proj3=proj3, quadwts=quadwts,
                dt=self.cfg.get('solver-time-integrator', 'dt')
            )

            # Apply the sensor to estimate the required artificial viscosity
            self.kernels['riemanndifference'] = lambda: self._be.kernel(
                'riemanndifference', tplargs=tplargs, dims=[self.neles],
                u=self.scal_upts_inb, plocu=self.ploc_at('upts'), usmats=self.ele_smat_at('upts'),
                uf=self._scal_fpts_cpy, fsmats=self.ele_smat_at('fpts'), divf=self._scal_upts_cpy
            )



    def generateDiffMat(self):
        p = self.basis.order
        M = np.zeros((p+1, p+2))

        solpts = self.basis.upts[:p+1,0] # Hack
        rdpts = np.zeros(p+2)
        rdpts[0] = -1.
        rdpts[-1] = 1.
        for i in range(p):
            rdpts[i+1] = 0.5*(solpts[i] + solpts[i+1])
        
        for i in range(p+2):
            vals = np.zeros(p+2)
            vals[i] = 1.
            dlag = interpolate.lagrange(rdpts, vals).deriv()
            for j in range(p+1):
                M[j,i] = dlag(solpts[j])
        return M

    def generateProjectionMats(self, pstages):
        # p1pts = np.array([-0.5773502, 0.5773502])
        # p2pts = np.array([-0.7745966, 0.0, 0.7745966])
        # p3pts = np.array([-0.8611363, -0.3399810, -0.3399810, -0.8611363])

        upts =  self.basis.upts
        solpts = self.basis.upts[:self.basis.order+1,0] # Hack
        p = self.basis.order
        if self.ndims == 2:
            # # Swap ordering of points for legvander2d
            # rsupts = np.reshape(upts, (self.basis.order+1, self.basis.order+1, 2))
            # rsupts = rsupts.swapaxes(0,1)
            # rsupts = np.reshape(rsupts, (-1, 2))

            # # Generate inverse VDM matrix for calculating Legendre modes
            # vdm  = np.polynomial.legendre.legvander2d(rsupts[:,0], rsupts[:,1], [self.basis.order]*2)
            # invvdm = np.linalg.inv(vdm)
            vdm = self.basis.ubasis.vdm.T
            invvdm = np.linalg.inv(vdm)
            
            # Generate VDM matrices
            vdm1 = copy.copy(vdm)
            vdm2 = copy.copy(vdm)
            vdm3 = copy.copy(vdm)

            modesfac1 = np.zeros((p+1, p+1))
            modesfac2 = np.zeros((p+1, p+1))
            modesfac3 = np.zeros((p+1, p+1))

            modesfac1[:pstages[0]+1, :pstages[0]+1] = 1.0
            modesfac2[:pstages[1]+1, :pstages[1]+1] = 1.0
            modesfac3[:pstages[2]+1, :pstages[2]+1] = 1.0

            modesfac1 = np.reshape(modesfac1, (-1))
            modesfac2 = np.reshape(modesfac2, (-1))
            modesfac3 = np.reshape(modesfac3, (-1))

            # Zero higher order modes
            for j in range(self.nupts):
                vdm1[:, j] *= modesfac1[j]
                vdm2[:, j] *= modesfac2[j]
                vdm3[:, j] *= modesfac3[j]

            # vdm1[:,(pstages[0]+1)**self.ndims:] *= 0
            # vdm2[:,(pstages[1]+1)**self.ndims:] *= 0
            # vdm3[:,(pstages[2]+1)**self.ndims:] *= 0

            # Create projection matrices
            proj1 = np.matmul(vdm1, invvdm)
            proj2 = np.matmul(vdm2, invvdm)
            proj3 = np.matmul(vdm3, invvdm)

            proj1[np.abs(proj1) < 1e-7] = 0.0
            proj2[np.abs(proj2) < 1e-7] = 0.0
            proj3[np.abs(proj3) < 1e-7] = 0.0

            # Create differentiation matrices
            nptsx = self.basis.order+1
            nptsy = self.basis.order+1
            npts = nptsx*nptsy

            Mxi = np.zeros((npts, npts))
            Meta = np.zeros((npts, npts))

            for i in range(nptsx):
                for j in range(nptsy):
                    vals = np.zeros(nptsx)
                    vals[i] = 1.
                    lag_xi = interpolate.lagrange(solpts, vals)
                    dlag_xi = lag_xi.deriv()
                    vals = np.zeros(nptsy)
                    vals[j] = 1.
                    lag_eta = interpolate.lagrange(solpts, vals)
                    dlag_eta = lag_eta.deriv()
                    for di in range(nptsx):
                        for dj in range(nptsy):
                            x = solpts[di]
                            y = solpts[dj]
                            Mxi[di+dj*nptsx, i+j*nptsx] = dlag_xi(x)*lag_eta(y)
                            Meta[di+dj*nptsx, i+j*nptsx] = lag_xi(x)*dlag_eta(y)

            # Get quadrature weights
            eletype = 'quad'
            rule = self.cfg.get('solver-elements-{0}'.format(eletype), 'soln-pts')
            wts = get_quadrule(eletype, rule=rule, npts=len(self.basis.upts), qdeg=self.basis.order, flags=None).wts

            return [Mxi, Meta, None, proj1, proj2, proj3, wts]
        elif self.ndims == 3:
            raise NotImplementedError()
            #invvdm = np.linalg.inv(np.polynomial.legendre.legvander3d(upts[:,0], upts[:,1], upts[:,2], [self.basis.order]*3))

        invvdm[np.abs(invvdm) < 1e-7] = 0.0

        return [vdm, invvdm] 
