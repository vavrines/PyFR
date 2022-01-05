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
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.calcminentropy')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs_inv = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float),
                       viscous=False)

        tplargs_vis = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float),
                       viscous=True)

        if 'flux' in self.antialias:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_inv, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_vis, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
        else:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_inv, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_vis, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )

        ubdegs = np.array([max(dd) for dd in self.basis.ubasis.degrees])
        ffac = np.exp(-ubdegs**2)
        dt = self.cfg.get('solver-time-integrator', 'dt')
        niters = int(self.cfg.get('solver', 'filter-iterations', 20))
        mean_mode_value = (self.basis.ubasis.invvdm.T @ np.ones_like(self.basis.upts[:,0]))[0]
        dtol = float(self.cfg.get('solver', 'filter-dtol', 1e-8))
        ptol = float(self.cfg.get('solver', 'filter-ptol', 1e-8))
        etol = float(self.cfg.get('solver', 'filter-etol', 1e-4))

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nupts=self.nupts, nfpts=self.nfpts,
                       c=self.cfg.items_as('constants', float), 
                       order=self.basis.order, ffac=ffac, dt=dt,
                       srcex=self._src_exprs, niters=niters,
                       dtol=dtol, ptol=ptol, etol=etol, 
                       mean_mode_value=mean_mode_value)

        plocupts = self.ploc_at('upts') if self._ploc_in_src_exprs else None

        self.kernels['negdivconf'] = lambda: self._be.kernel(
            'negdivconffilter', tplargs=tplargs,
            dims=[self.neles], tdivtconf_vis=self._scal_upts_cpy, tdivtconf_inv=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=self.scal_upts_inb,
            entmin=self.entmin, vdm=self.vdm, invvdm=self.invvdm
        )

        self.kernels['calcentropy'] = lambda: self._be.kernel(
            'calcentropy', tplargs=tplargs, dims=[self.neles],
            u=self.scal_upts_inb, entmin=self.entmin
        )

        self.kernels['calcminentropy'] = lambda: self._be.kernel(
            'calcminentropy', tplargs=tplargs, dims=[self.neles],
            entmin=self.entmin, entmin_int=self.entmin_int
        )




