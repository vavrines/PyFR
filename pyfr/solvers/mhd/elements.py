# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.solvers.euler.elements import BaseFluidElements
import numpy as np

class MHDElements(BaseAdvectionElements):
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

    @staticmethod
    def validate_formulation(ctrl):
        shock_capturing = ctrl.cfg.get('solver', 'shock-capturing', 'none')
        if ctrl.formulation == 'dual' and shock_capturing == 'entropy-filter':
            raise ValueError('Entropy filtering not compatible with '
                             'dual time stepping.')

        ctrlvardt = ctrl.controller_has_variable_dt
        if ctrlvardt and shock_capturing == 'entropy-filter':
            raise ValueError('Entropy filtering not compatible with '
                             'adaptive time stepping.')

    @property
    def _soln_in_src_exprs(self):
        return True

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mhd.kernels.negdivconfmhd')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.powellsource')

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs
        solnsrc = True

        # Transformed to physical divergence kernel + source term
        plocupts = self.ploc_at('upts') if plocsrc else None
        solnupts = self._scal_upts_cpy if solnsrc else None

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'srcex': self._src_exprs
        }

        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconfmhd', tplargs=srctplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts
        )

        self.kernels['powellsource'] = lambda fout: self._be.kernel(
            'powellsource', tplargs=srctplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts
        )

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tfluxlin')


        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
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
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
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
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )


        # Can elide shock-capturing at p = 0
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')

        # Modified entropy filtering method using specific physical
        # entropy (without operator splitting for Navier-Stokes)
        # doi:10.1016/j.jcp.2022.111501
        if shock_capturing == 'entropy-filter' and self.basis.order != 0:
            self._be.pointwise.register(
                'pyfr.solvers.mhd.kernels.entropylocal'
            )
            self._be.pointwise.register(
                'pyfr.solvers.mhd.kernels.entropyfilter'
            )

            # Template arguments
            eftplargs = {
                'ndims': self.ndims,
                'nupts': self.nupts,
                'nfpts': self.nfpts,
                'nvars': self.nvars,
                'nfaces': self.nfaces,
                'c': self.cfg.items_as('constants', float),
                'order': self.basis.order
            }

            # Check to see if running collocated solution/flux points
            m0 = self.basis.m0
            mrowsum = np.max(np.abs(np.sum(m0, axis=1) - 1.0))
            if np.min(m0) < -1e-8 or mrowsum > 1e-8:
                raise ValueError('Entropy filter requires flux points to be a '
                                 'subset of solution points or a convex '
                                 'combination thereof.')

            # Minimum density/pressure constraints
            eftplargs['d_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'd-min', 1e-6)
            eftplargs['p_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'p-min', 1e-6)

            # Entropy tolerance
            eftplargs['e_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'e-tol', 1e-6)

            # Hidden kernel parameters
            eftplargs['f_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'f-tol', 1e-4)
            eftplargs['ill_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                     'ill-tol', 1e-6)
            eftplargs['niters'] = self.cfg.getfloat('solver-entropy-filter',
                                                    'niters', 20)

            # Precompute basis orders for filter
            ubdegs = self.basis.ubasis.degrees
            eftplargs['ubdegs'] = [int(max(dd)) for dd in ubdegs]
            eftplargs['order'] = self.basis.order

            # Compute local entropy bounds
            self.kernels['local_entropy'] = lambda uin: self._be.kernel(
                'entropylocal', tplargs=eftplargs, dims=[self.neles],
                u=self.scal_upts[uin], entmin_int=self.entmin_int
            )

            # Apply entropy filter
            self.kernels['entropy_filter'] = lambda uin: self._be.kernel(
                'entropyfilter', tplargs=eftplargs, dims=[self.neles],
                u=self.scal_upts[uin], entmin_int=self.entmin_int,
                vdm=self.vdm, invvdm=self.invvdm
            )