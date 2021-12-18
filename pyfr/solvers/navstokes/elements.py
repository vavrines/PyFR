# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements
import numpy as np

class NavierStokesElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.negdivconffilter')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.calcentropy')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float))

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
            entmin=self.entmin
        )

        # Apply the sensor to estimate the required artificial viscosity
        self.kernels['calcentropy'] = lambda: self._be.kernel(
            'calcentropy', tplargs=tplargs, dims=[self.neles],
            u=self.scal_upts_inb, entmin=self.entmin, entmin_cpy=self.entmin_cpy
        )


