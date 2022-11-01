# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

import numpy as np
from math import gamma as gamma_func

def reflect2D(side, Nr, Nt, Nz=1):
      N = Nr*Nt*Nz
      pairs = np.zeros(N, dtype=int)

      assert Nt % 4 == 0, 'Nt must be div4 for wall BCs.'

      if side == 'lr':
            tpairs = np.linspace(0, Nt-1, Nt, dtype=int)[::-1]
            tpairs = np.roll(tpairs, Nt//2 + 1) 
      elif side == 'ud':
            tpairs = np.linspace(0, Nt-1, Nt, dtype=int)[::-1]
            tpairs = np.roll(tpairs, 1)
      elif side == 'diagr':
            tpairs = np.linspace(0, Nt-1, Nt, dtype=int)[::-1]
            tpairs = np.roll(tpairs, Nt//4 + 1)

      n = 0
      for i in range(Nz):
            for j in range(Nr):
                  for k in range(Nt):
                        pairs[n + k] = n + tpairs[k]
                  n += Nt
      return pairs

class BGKIntInters(BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.intcflux')


        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat
        )


class BGKMPIInters(BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.mpicflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat
        )


class BGKBaseBCInters(BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.bccflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.niters = self.cfg.getint('solver', 'niters')
        delta = self.cfg.getint('solver', 'delta')
        lam = 1.0/gamma_func(delta/2.0) if delta else 0.0

                # Get reflections for wall BCs
        Nr = self.cfg.getint('solver', 'Nr')
        Nt = self.cfg.getint('solver', 'Nt')
        Nz = self.cfg.getint('solver', 'Nz') if delta else 1

        self.LRidxs = reflect2D('lr', Nr, Nt, Nz)
        self.UDidxs = reflect2D('ud', Nr, Nt, Nz)
        self.DRidxs = reflect2D('diagr', Nr, Nt, Nz)

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       c=self.c, u=self.u, bctype=self.type, niters=self.niters,
                       pi=np.pi, delta=delta, lam=lam, LRidxs=self.LRidxs, 
                       UDidxs=self.UDidxs, DRidxs=self.DRidxs)
        
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

class BGKDiffuseBCInters(BGKBaseBCInters):
    type = 'diffuse'

class BGKSpecularBCInters(BGKBaseBCInters):
    type = 'specular'


