# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements
import numpy as np

class MHDElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'
    formulations = ['std', 'dual']
    privarmap = {2: ['rho', 'u', 'v', 'Bx', 'By', 'divB', 'p'],
                 3: ['rho', 'u', 'v', 'w', 'Bx', 'By', 'Bz', 'divB', 'p']}

    convarmap = {2: ['rho', 'rhou', 'rhov', 'Bx', 'By', 'divB', 'E'],
                 3: ['rho', 'rhou', 'rhov', 'rhow', 'Bx', 'By', 'Bz', 'divB', 'E']}

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('B', ['Bx', 'By']),
            ('pressure', ['p']),
            ('divB', ['divB'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('B', ['Bx', 'By', 'Bz']),
            ('pressure', ['p']),
            ('divB', ['divB'])]
    }

    @staticmethod
    def pri_to_con(pris, cfg):
        if len(pris) == 7:
            ndims = 2
        elif len(pris) == 9:
            ndims = 3
        else:
            raise ValueError('Unknown length of vars.', pris)
        # Density, pressure, velocity field, and magnetic field
        rho, divb, p = pris[0], pris[-2], pris[-1]
        vf, Bf = list(pris[1:ndims+1]), list(pris[ndims+1:2*ndims+1])

        # Multiply velocity components by rho
        rhovf = [rho*v for v in vf]

        # Squared velocity and magnetic fields
        vf2 = sum(v*v for v in vf)
        bf2 = sum(b*b for b in Bf)

        # Compute the energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*vf2 + 0.5*bf2

        return [rho] + rhovf + Bf + [divb] + [E]

    @staticmethod
    def con_to_pri(cons, cfg):
        if len(cons) == 7:
            ndims = 2
        elif len(cons) == 9:
            ndims = 3
        else:
            raise ValueError('Unknown length of vars.', cons)
        # Density, energy, momentum field, and magnetic field
        rho, divb, E = cons[0], cons[-2], cons[-1]
        rhovf, Bf = list(cons[1:ndims+1]), list(cons[ndims+1:2*ndims+1])

        # Divide momentum components by rho
        vf = [rhov/rho for rhov in rhovf]

        # Squared velocity and magnetic fields
        vf2 = sum(v*v for v in vf)
        bf2 = sum(b*b for b in Bf)

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*vf2 - 0.5*bf2)

        return [rho] + vf + Bf + [divb] + [p]

    @property
    def _soln_in_src_exprs(self):
        return True

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.negdivconffilter')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.calcentropy')


        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float),
                       srcex=self._src_exprs)

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
        else:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )


        ubdegs = np.array([max(dd) for dd in self.basis.ubasis.degrees])
        ffac = np.exp(-ubdegs**2)
        dt = self.cfg.get('solver-time-integrator', 'dt')
        niters = int(self.cfg.get('solver', 'filter-iterations', 10))
        alpha = float(self.cfg.get('solver', 'filter-alpha', 1.0))
        mean_mode_value = (self.basis.ubasis.invvdm.T @ np.ones_like(self.basis.upts[:,0]))[0]
        dtol = 1e-8
        ptol = 1e-8
        etol = 1e-8

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nupts=self.nupts,
                       c=self.cfg.items_as('constants', float), 
                       order=self.basis.order, ffac=ffac, dt=dt,
                       vdm=self.basis.ubasis.vdm.T,
                       invvdm=self.basis.ubasis.invvdm.T,
                       srcex=self._src_exprs, niters=niters,
                       dtol=dtol, ptol=ptol, etol=etol, 
                       alpha=alpha, mean_mode_value=mean_mode_value)

        plocupts = self.ploc_at('upts') if self._ploc_in_src_exprs else None

        self.kernels['negdivconf'] = lambda: self._be.kernel(
            'negdivconffilter', tplargs=tplargs,
            dims=[self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=self.scal_upts_inb,
            ent_min=self.artvisc
        )

        # Apply the sensor to estimate the required artificial viscosity
        self.kernels['calcentropy'] = lambda: self._be.kernel(
            'calcentropy', tplargs=tplargs, dims=[self.neles],
            u=self.scal_upts_inb, ent_min=self.artvisc
        )
