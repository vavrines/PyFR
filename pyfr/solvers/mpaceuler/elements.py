# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements


class BaseMPACFluidElements:
    privarmap = {2: ['p', 'u', 'v', 'phi'],
                 3: ['p', 'u', 'v', 'w', 'phi']}

    pasvarmap = {2: ['rho'], 3: ['rho']}

    convarmap = {2: ['p', 'rhou', 'rhov', 'phi'],
                 3: ['p', 'rhou', 'rhov', 'rhow', 'phi']}

    dualcoeffs = {2: ['rhou', 'rhov'],
                  3: ['rhou', 'rhov', 'rhow']}

    visvarmap = {
        2: [('velocity', ['u', 'v']),
            ('pressure', ['p']),
            ('phase', ['phi'])],
        3: [('velocity', ['u', 'v', 'w']),
            ('pressure', ['p']),
            ('phase', ['phi'])]
    }

    @staticmethod
    def pri_to_con(pris, cfg, passive):
        rho = passive[0]
        p, phi = pris[0], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-1]]
        return [p] + rhovs + [phi]

    @staticmethod
    def con_to_pri(cons, cfg, passive):
        rho = passive[0]
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
                u=s(self.scal_upts[uin], c), q=s(self._pass_upts, c), 
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
                u=s(self.scal_upts[uin], l), q=s(self._pass_upts, l),
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
            u=self.scal_upts[uin], p=self._pass_upts
        )

        # Interpolation from elemental points
        self.kernels['disq'] = lambda uin: self._be.kernel(
            'mul', self.opmat('M0'), self._pass_upts[uin],
            out=self._pass_fpts
        )