# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

import numpy as np

class BGKIntInters(BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.intcflux')


        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c, u=self.u)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class BGKMPIInters(BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.mpicflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c, u=self.u)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class BGKBaseBCInters(BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.bccflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.niters = self.cfg.getint('solver', 'niters')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c, u=self.u, bctype=self.type, niters=self.niters,
                       pi=np.pi)
        
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, fl=self._scal_lhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat, M=self.M, **self._external_vals
        )

class BGKFreeBCInters(BGKBaseBCInters):
    type = 'free'
    cflux_state = 'ghost'

class BGKFixedBCInters(BGKBaseBCInters):
    type = 'fixed'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )

class BGKSlipWallBCInters(BGKBaseBCInters):
    type = 'slip-wall'

class BGKNoSlipWallBCInters(BGKBaseBCInters):
    type = 'no-slip-wall'


