# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np
from scipy import interpolate

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
    def pri_to_con(pris, cfg):
        rho, p = pris[0], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-1]]

        # Compute the energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[1:-1])

        return [rho] + rhovs + [E]

    @staticmethod
    def con_to_pri(cons, cfg):
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

        if shock_capturing == 'riemann-difference':
            # Register the kernels
            self._be.pointwise.register('pyfr.solvers.euler.kernels.rdshocksensor')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.riemanndifference')

            self.kernels['copy_soln_at_fpts'] = lambda: self._be.kernel(
                    'copy', self._scal_fpts_cpy, self._scal_fpts
                )

            # Sense shocks using density convergence
            shockvar = 0

            # Template arguments
            sftplargs = dict(
                nvars=self.nvars, ndims=self.ndims, nupts=self.nupts, svar=shockvar,
                c=self.cfg.items_as('solver-riemann-difference', float),
                order=self.basis.order  
            )

            # Shockcell is an elementwise value (0,1) dictating if there is a shock in the cell (1)
            self.shockcell = self._be.matrix((1, self.neles), tags={'align'}, initval=np.zeros((1, self.neles)))

            # Currently the sensor is always on and applies filter everywhere
            self.kernels['rdshocksensor'] = lambda: self._be.kernel(
               'rdshocksensor', tplargs=sftplargs, dims=[self.neles],
               u=self.scal_upts_inb, shockcell=self.shockcell
            )

            diffmat = self.generateDiffMat()
            rsolver = self.cfg.get('solver-interfaces', 'riemann-solver') + '3d'

            # Template arguments
            tplargs = dict(
                nvars=self.nvars, nupts=self.nupts, nfpts=self.nfpts, ndims=self.ndims,
                c=self.cfg.items_as('constants', float),order=self.basis.order, 
                diffmat=diffmat, rsolver=rsolver
            )

            # Apply the sensor to estimate the required artificial viscosity
            self.kernels['riemanndifference'] = lambda: self._be.kernel(
                'riemanndifference', tplargs=tplargs, dims=[self.neles],
                u=self.scal_upts_inb, plocu=self.ploc_at('upts'), urcpdjac=self.rcpdjac_at('upts'), usmats=self.ele_smat_at('upts'),
                uf=self._scal_fpts_cpy, frcpdjac=self.rcpdjac_at('fpts'), fsmats=self.ele_smat_at('fpts'),
                divf=self.scal_upts_outb
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
