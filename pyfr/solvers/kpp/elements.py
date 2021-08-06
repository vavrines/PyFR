# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements


class KPPElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'
    formulations = ['std', 'dual']
    privarmap = {2: ['u'],
                 3: ['u']}

    convarmap = {2: ['u'],
                 3: ['u']}

    visvarmap = {
        2: [('u', ['u'])],
        3: [('u', ['u'])]
    }

    @staticmethod
    def pri_to_con(pris, cfg):
        return pris

    @staticmethod
    def con_to_pri(cons, cfg):
        return cons

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.kpp.kernels.tflux_inv')
        self._be.pointwise.register('pyfr.solvers.kpp.kernels.tflux_vis')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.get_du')
        self._be.pointwise.register('pyfr.solvers.baseadvec.kernels.negdivconf_f')
        self._be.pointwise.register('pyfr.solvers.baseadvec.kernels.negdivconf_b')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        dt_rev = float(self.cfg.get('solver-artificial-viscosity', 'dt_rev', self.cfg.get('solver-time-integrator', 'dt')))
        mu_max = float(self.cfg.get('solver-artificial-viscosity', 'mu_max'))
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float),
                       srcex=self._src_exprs, dt_rev=dt_rev,
                       mu_max=mu_max)

        if 'flux' in self.antialias:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux_inv', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux_vis', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
        else:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux_inv', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux_vis', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )



        self.kernels['negdivconf_f'] = lambda: self._be.kernel(
            'negdivconf_f', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb
        )
        self.kernels['negdivconf_b'] = lambda: self._be.kernel(
            'negdivconf_b', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb
        )

        self.kernels['get_du'] = lambda: self._be.kernel(
            'get_du', tplargs=tplargs,
            dims=[self.nupts, self.neles], u_orig=self._scal_upts_cpy, 
            u_rev=self.scal_upts_inb
        )


