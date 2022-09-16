# -*- coding: utf-8 -*-

import numpy as np

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.mpeuler.elements import BaseMPFluidElements


class MPNavierStokesElements(BaseMPFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    @property
    def _scratch_bufs(self):

        bufs = {'scal_fpts', 'vect_fpts', 'vect_upts', 'vect_upts_cpy',
                'scal_upts_cpy'}

        if 'flux' in self.antialias:
            bufs |= {'vect_qpts_cpy', 'scal_qpts', 'vect_qpts'}

        return bufs

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
        self._be.pointwise.register('pyfr.solvers.mpnavstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mpnavstokes.kernels.tfluxlin')
        self._be.pointwise.register('pyfr.solvers.mpnavstokes.kernels.negdivconf', force=True)

        # Handle shock capturing and Sutherland's law
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nspec': self.nspec,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'shock_capturing': shock_capturing,
            'visc_corr': visc_corr
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        av = self.artvisc

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), grad=s(self._vect_upts_cpy, l),
                f=s(self._vect_upts, c), artvisc=s(av, c), 
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), grad=s(self._vect_qpts_cpy, l),
                f=s(self._vect_qpts, c), artvisc=s(av, c), 
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), grad=s(self._vect_upts_cpy, l),
                f=s(self._vect_upts, l), artvisc=s(av, l), 
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), grad=s(self._vect_qpts_cpy, l),
                f=s(self._vect_qpts, l), artvisc=s(av, l), 
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )

        self.kernels['copy_grad'] = lambda: self._be.kernel(
            'copy', self._vect_upts_cpy, self._vect_upts
        )

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nspec': self.nspec,
            'srcex': self._src_exprs
        }

        # Transformed to physical divergence kernel + source term
        plocupts = self.ploc_at('upts') if plocsrc else None

        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconf', tplargs=srctplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, 
            u=self._scal_upts_cpy, grad=self._vect_upts_cpy,
        )

        # Re-register kernels for gradient calcs
        kernel = self._be.kernel
        slicem = self._slice_mat

        # Mesh regions
        regions = self._mesh_regions

        if self.basis.order > 0:
            self.kernels['tgradpcoru_upts'] = lambda uin: kernel(
                'mul', self.opmat('M4 - M6*M0'), self.scal_upts[uin],
                out=self._vect_upts_cpy
            )
        self.kernels['tgradcoru_upts'] = lambda: kernel(
            'mul', self.opmat('M6'), self._vect_fpts.slice(0, self.nfpts),
            out=self._vect_upts_cpy, beta=float(self.basis.order > 0)
        )

        if 'curved' in regions:
            self.kernels['gradcoru_upts_curved'] = lambda: kernel(
                'gradcoru', tplargs=tplargs,
                dims=[self.nupts, regions['curved']],
                gradu=slicem(self._vect_upts_cpy, 'curved'),
                smats=self.curved_smat_at('upts'),
                rcpdjac=self.rcpdjac_at('upts', 'curved')
            )

        if 'linear' in regions:
            self.kernels['gradcoru_upts_linear'] = lambda: kernel(
                'gradcorulin', tplargs=tplargs,
                dims=[self.nupts, regions['linear']],
                gradu=slicem(self._vect_upts_cpy, 'linear'),
                upts=self.upts, verts=self.ploc_at('linspts', 'linear')
            )

        def gradcoru_fpts():
            nupts, nfpts = self.nupts, self.nfpts
            vupts, vfpts = self._vect_upts_cpy, self._vect_fpts

            # Exploit the block-diagonal form of the operator
            muls = [kernel('mul', self.opmat('M0'),
                           vupts.slice(i*nupts, (i + 1)*nupts),
                           vfpts.slice(i*nfpts, (i + 1)*nfpts))
                    for i in range(self.ndims)]

            return self._be.unordered_meta_kernel(muls)

        self.kernels['gradcoru_fpts'] = gradcoru_fpts

        if 'flux' in self.antialias and self.basis.order > 0:
            def gradcoru_qpts():
                nupts, nqpts = self.nupts, self.nqpts
                vupts, vqpts = self._vect_upts_cpy, self._vect_qpts_cpy

                # Exploit the block-diagonal form of the operator
                muls = [self._be.kernel('mul', self.opmat('M7'),
                                        vupts.slice(i*nupts, (i + 1)*nupts),
                                        vqpts.slice(i*nqpts, (i + 1)*nqpts))
                        for i in range(self.ndims)]

                return self._be.unordered_meta_kernel(muls)

            self.kernels['gradcoru_qpts'] = gradcoru_qpts
