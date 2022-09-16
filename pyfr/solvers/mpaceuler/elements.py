# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
from symbol import classdef


class BaseMPACFluidElements:
    @classmethod
    def privarmap(cls, cfg=None):
        return {2: ['p', 'u', 'v', 'phi'],
                3: ['p', 'u', 'v', 'w', 'phi']}

    @classmethod
    def pasvarmap(cls, cfg=None):
        return {2: ['rho'], 3: ['rho']}

    @classmethod
    def convarmap(cls, cfg=None):
        return {2: ['p', 'rhou', 'rhov', 'phi'],
                3: ['p', 'rhou', 'rhov', 'rhow', 'phi']}

    @classmethod
    def dualcoeffs(cls, cfg=None):
        return {2: ['rhou', 'rhov', 'phi'],
                3: ['rhou', 'rhov', 'rhow', 'phi']}

    @classmethod
    def visvarmap(cls, cfg=None):
        return {2: [('velocity', ['u', 'v']),
                    ('pressure', ['p']),
                    ('phase', ['phi'])],
                3: [('velocity', ['u', 'v', 'w']),
                    ('pressure', ['p']),
                    ('phase', ['phi'])]
               }

    @property
    def _scratch_bufs(self):
        if 'flux' in self.antialias:
            bufs = {'scal_fpts', 'scal_qpts', 'vect_qpts'}
        else:
            bufs = {'scal_fpts', 'vect_upts'}

        if self._soln_in_src_exprs:
            bufs |= {'scal_upts_cpy'}

        bufs |= {'pasv_upts', 'pasv_fpts'}

        return bufs

    @staticmethod
    def pri_to_con(pris, cfg, pasv):
        rho = pasv[0]
        p, phi = pris[0], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-1]]
        return [p] + rhovs + [phi]

    @staticmethod
    def con_to_pri(cons, cfg, pasv):
        rho = pasv[0]
        p, phi = cons[0], cons[-1]

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[1:-1]]

        return [p] + vs + [phi]

    @staticmethod
    def validate_formulation(controller):
        if controller.formulation != 'dual':
            raise ValueError('System not compatible with time stepping '
                             'formulation.')


class MPACEulerElements(BaseMPACFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mpaceuler.kernels.density')
        self._be.pointwise.register('pyfr.solvers.mpaceuler.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mpaceuler.kernels.tfluxlin')
        self._be.pointwise.register(
            'pyfr.solvers.mpaceuler.kernels.negdivconf', force=True
        )

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'npass': self.npass,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), q=s(self.pasv_upts, c), 
                f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), q=s(self.pasv_upts, l),
                f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )

        self.kernels['density'] = lambda uin: self._be.kernel(
            'density', tplargs=tplargs, dims=[self.nupts, self.neles],
            u=self.scal_upts[uin], q=self.pasv_upts
        )

        # Interpolation from elemental points
        self.kernels['disq'] = lambda uin: self._be.kernel(
            'mul', self.opmat('M0'), self.pasv_upts,
            out=self._pasv_fpts
        )

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs
        solnsrc = self._soln_in_src_exprs
        pasvsrc = self._pasv_in_src_exprs

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'npass': self.npass,
            'srcex': self._src_exprs
        }

        # Transformed to physical divergence kernel + source term
        plocupts = self.ploc_at('upts') if plocsrc else None
        solnupts = self._scal_upts_cpy if solnsrc else None
        pasvupts = self.pasv_upts if pasvsrc else None
        
        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconf', tplargs=srctplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts, 
            q=pasvupts,
        )
