# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

class FluidIntIntersMixin:
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        self._be.pointwise.register('pyfr.solvers.euler.kernels.intcent')

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'intcent', tplargs={}, dims=[self.ninters],
                entmin_lhs=self._entmin_lhs, entmin_rhs=self._entmin_rhs
            )


class FluidMPIIntersMixin:
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        self._be.pointwise.register('pyfr.solvers.euler.kernels.mpicent')

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'mpicent', tplargs={}, dims=[self.ninters],
                entmin_lhs=self._entmin_lhs, entmin_rhs=self._entmin_rhs
            )


class EulerIntInters(FluidIntIntersMixin, BaseAdvectionIntInters):
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        self._be.pointwise.register('pyfr.solvers.euler.kernels.intcflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class EulerMPIInters(FluidMPIIntersMixin, BaseAdvectionMPIInters):
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        self._be.pointwise.register('pyfr.solvers.euler.kernels.mpicflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class EulerBaseBCInters(BaseAdvectionBCInters):
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self._be.pointwise.register('pyfr.solvers.euler.kernels.bccflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c, bctype=self.type, ninters=self.ninters)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self._be.pointwise.register('pyfr.solvers.euler.kernels.bccent')

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )


class EulerSupInflowBCInters(EulerBaseBCInters):
    type = 'sup-in-fa'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class EulerSupOutflowBCInters(EulerBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class EulerCharRiemInvBCInters(EulerBaseBCInters):
    type = 'char-riem-inv'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class EulerSlpAdiaWallBCInters(EulerBaseBCInters):
    type = 'slp-adia-wall'
