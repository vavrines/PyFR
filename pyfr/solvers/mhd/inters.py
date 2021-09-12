# -*- coding: utf-8 -*-

import numpy as np

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)


class TplargsMixin(object):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = self.cfg.get('solver', 'viscosity-correction')
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        mu_max = float(self.cfg.get('solver-artificial-viscosity', 'mu_max'))
        lambda_min = self.cfg.get('solver-interfaces', 'lambda-min', 0.0)
        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver, visc_corr=visc_corr,
                             shock_capturing=shock_capturing, c=self.c,
                             mu_max=mu_max, lambda_min=lambda_min)


class MHDIntInters(TplargsMixin, BaseAdvectionDiffusionIntInters):
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.mhd.kernels.intconu')
        be.pointwise.register('pyfr.solvers.mhd.kernels.intcflux_vis')
        be.pointwise.register('pyfr.solvers.mhd.kernels.intcflux_inv_f')
        be.pointwise.register('pyfr.solvers.mhd.kernels.intcflux_inv_b')

        if abs(self.c['ldg-beta']) == 0.5:
            self.kernels['copy_fpts'] = lambda: ComputeMetaKernel(
                [ele.kernels['_copy_fpts']() for ele in elemap.values()]
            )

        self.kernels['con_u'] = lambda: be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs,
            ulout=self._vect_lhs, urout=self._vect_rhs
        )
        self.kernels['comm_flux_vis'] = lambda: be.kernel(
            'intcflux_vis', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )
        self.kernels['comm_flux_inv_f'] = lambda: be.kernel(
            'intcflux_inv_f', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )
        self.kernels['comm_flux_inv_b'] = lambda: be.kernel(
            'intcflux_inv_b', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class MHDMPIInters(TplargsMixin, BaseAdvectionDiffusionMPIInters):
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.mhd.kernels.mpiconu')
        be.pointwise.register('pyfr.solvers.mhd.kernels.mpicflux_vis')
        be.pointwise.register('pyfr.solvers.mhd.kernels.mpicflux_inv_f')
        be.pointwise.register('pyfr.solvers.mhd.kernels.mpicflux_inv_b')

        self.kernels['con_u'] = lambda: be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs, ulout=self._vect_lhs
        )
        self.kernels['comm_flux_vis'] = lambda: be.kernel(
            'mpicflux_vis', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )
        self.kernels['comm_flux_inv_f'] = lambda: be.kernel(
            'mpicflux_inv_f', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )
        self.kernels['comm_flux_inv_b'] = lambda: be.kernel(
            'mpicflux_inv_b', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class MHDBaseBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type
        self._tplargs['bccfluxstate'] = self.cflux_state

        be.pointwise.register('pyfr.solvers.mhd.kernels.bcconu')
        be.pointwise.register('pyfr.solvers.mhd.kernels.bccflux_vis')
        be.pointwise.register('pyfr.solvers.mhd.kernels.bccflux_inv_f')
        be.pointwise.register('pyfr.solvers.mhd.kernels.bccflux_inv_b')

        self.kernels['con_u'] = lambda: be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._vect_lhs, nlin=self._norm_pnorm_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux_vis'] = lambda: be.kernel(
            'bccflux_vis', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs,
            gradul=self._vect_lhs, magnl=self._mag_pnorm_lhs,
            nl=self._norm_pnorm_lhs, artviscl=self._artvisc_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux_inv_f'] = lambda: be.kernel(
            'bccflux_inv_f', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs,
            gradul=self._vect_lhs, magnl=self._mag_pnorm_lhs,
            nl=self._norm_pnorm_lhs, artviscl=self._artvisc_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux_inv_b'] = lambda: be.kernel(
            'bccflux_inv_b', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs,
            gradul=self._vect_lhs, magnl=self._mag_pnorm_lhs,
            nl=self._norm_pnorm_lhs, artviscl=self._artvisc_lhs,
            **self._external_vals
        )


class MHDSlpAdiaWallBCInters(MHDBaseBCInters):
    type = 'diode'
    cflux_state = None


class MHDFreeBCInters(MHDBaseBCInters):
    type = 'free'
    cflux_state = 'ghost'


class MHDFixedBCInters(MHDBaseBCInters):
    type = 'fixed'
    cflux_state = None

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        if self.ndims == 2:
            tplc = self._exp_opts(
                ['rho', 'p', 'u', 'v', 'Bx', 'By'], lhs
            )
        else:
            tplc = self._exp_opts(
                ['rho', 'p', 'u', 'v', 'Bx', 'By', 'w', 'Bz'], lhs
            )
        self.c.update(tplc)