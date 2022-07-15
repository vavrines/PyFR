# -*- coding: utf-8 -*-

import numpy as np

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements

import numpy as np
import math

class NavierStokesElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        rho, *rhouvw = cons[:-1]
        grad_rho, *grad_rhouvw, grad_E = grad_cons

        # Divide momentum components by ρ
        uvw = [rhov / rho for rhov in rhouvw]

        # Velocity gradients: ∇u⃗ = 1/ρ·[∇(ρu⃗) - u⃗ ⊗ ∇ρ]
        grad_uvw = [(grad_rhov - v*grad_rho) / rho
                    for grad_rhov, v in zip(grad_rhouvw, uvw)]

        # Pressure gradient: ∇p = (γ - 1)·[∇E - 1/2*(u⃗·∇(ρu⃗) - ρu⃗·∇u⃗)]
        gamma = cfg.getfloat('constants', 'gamma')
        grad_p = grad_E - 0.5*(np.einsum('ijk,iljk->ljk', uvw, grad_rhouvw) +
                               np.einsum('ijk,iljk->ljk', rhouvw, grad_uvw))
        grad_p *= (gamma - 1)

        return [grad_rho] + grad_uvw + [grad_p]

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.tfluxlin')

        # Handle shock capturing and Sutherland's law
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        # Template parameters for the flux kernels
        tplargs_inv = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'shock_capturing': shock_capturing,
            'visc_corr': visc_corr
        }
        tplargs_vis = dict(tplargs_inv)
        tplargs_vis['viscous'] = True

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        av = self.artvisc

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved_inv'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs_inv, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('upts')
            )
            self.kernels['tdisf_curved_vis'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs_vis, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved_inv'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_inv, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('qpts')
            )
            self.kernels['tdisf_curved_vis'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs_vis, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear_inv'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs_inv, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.upts
            )
            self.kernels['tdisf_linear_vis'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs_vis, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear_inv'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs_inv, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.qpts
            )
            self.kernels['tdisf_linear_vis'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs_vis, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.qpts
            )

            
        if self.cfg.get('solver', 'shock-capturing', 'none') == 'entropy-filter':
            self._be.pointwise.register('pyfr.solvers.euler.kernels.entropylocal')
            self._be.pointwise.register('pyfr.solvers.euler.kernels.entropyfilter')

            # Entropy filtering not compatible with dual time
            self.formulations = ['std']

            # Minimum density/pressure constraints
            d_min = self.cfg.getfloat('solver-entropy-filter', 'd_min', 1e-6)
            p_min = self.cfg.getfloat('solver-entropy-filter', 'p_min', 1e-6)
            # Entropy tolerance
            e_tol = self.cfg.getfloat('solver-entropy-filter', 'e_tol', 1e-4)
            # Maximum filter strength (based on machine precision)
            eps = np.finfo(self._be.fpdtype).eps
            zeta_max = -math.log(eps)

            # Precompute basis orders for filter
            ubdegs2 = [max(dd)**2 for dd in self.basis.ubasis.degrees]

            eftplargs_inv = {
                'ndims': self.ndims, 'nupts': self.nupts, 'nfpts': self.nfpts,
                'nvars': self.nvars, 'c': self.cfg.items_as('constants', float),
                'd_min': d_min, 'p_min': p_min, 'e_tol': e_tol, 'zeta_max': zeta_max,
                'ubdegs2': ubdegs2, 'viscous': False
            }
            eftplargs_vis = dict(eftplargs_inv)
            eftplargs_vis['viscous'] = True

            # Compute local entropy bounds
            self.kernels['local_entropy'] = lambda uin: self._be.kernel(
                'entropylocal', tplargs=eftplargs_inv, dims=[self.neles],
                u=self.scal_upts[uin], entmin=self.entmin, 
                entmin_int=self.entmin_int
            )

            # Apply entropy filter
            self.kernels['filter_solution_inv'] = lambda uin: self._be.kernel(
                'entropyfilter', tplargs=eftplargs_inv, dims=[self.neles],
                u=self.scal_upts[uin], entmin=self.entmin,
                vdm=self.vdm, invvdm=self.invvdm
            )

            # Apply entropy filter
            self.kernels['filter_solution_vis'] = lambda uin: self._be.kernel(
                'entropyfilter', tplargs=eftplargs_vis, dims=[self.neles],
                u=self.scal_upts[uin], entmin=self.entmin,
                vdm=self.vdm, invvdm=self.invvdm
            )
