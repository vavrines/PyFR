# -*- coding: utf-8 -*-

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.quadrules import get_quadrule
import numpy as np

class BaseAdvectionDiffusionElements(BaseAdvectionElements):
    @property
    def _scratch_bufs(self):
        bufs = {'scal_fpts', 'vect_fpts', 'vect_upts'}

        if 'flux' in self.antialias:
            bufs |= {'scal_qpts', 'vect_qpts'}

        bufs |= {'scal_upts_cpy'}

        return bufs

    def set_backend(self, backend, nscalupts, nonce):
        super().set_backend(backend, nscalupts, nonce)

        kernel = self._be.kernel
        kernels = self.kernels

        # Register pointwise kernels
        self._be.pointwise.register(
            'pyfr.solvers.baseadvecdiff.kernels.gradcoru'
        )

        kernels['copy_soln'] = lambda: self._be.kernel(
                'copy', self._scal_upts_cpy, self.scal_upts_inb
            )
        
        kernels['_copy_fpts'] = lambda: kernel(
            'copy', self._vect_fpts.slice(0, self.nfpts), self._scal_fpts
        )
        kernels['tgradpcoru_upts'] = lambda: kernel(
            'mul', self.opmat('M4 - M6*M0'), self.scal_upts_inb,
            out=self._vect_upts
        )

        kernels['tgradcoru_upts'] = lambda: kernel(
            'mul', self.opmat('M6'), self._vect_fpts.slice(0, self.nfpts),
            out=self._vect_upts, beta=1.0
        )
        kernels['gradcoru_upts'] = lambda: kernel(
            'gradcoru', tplargs=dict(ndims=self.ndims, nvars=self.nvars),
            dims=[self.nupts, self.neles], smats=self.smat_at('upts'),
            rcpdjac=self.rcpdjac_at('upts'), gradu=self._vect_upts
        )

        def gradcoru_fpts():
            nupts, nfpts = self.nupts, self.nfpts
            vupts, vfpts = self._vect_upts, self._vect_fpts

            # Exploit the block-diagonal form of the operator
            muls = [kernel('mul', self.opmat('M0'),
                           vupts.slice(i*nupts, (i + 1)*nupts),
                           vfpts.slice(i*nfpts, (i + 1)*nfpts))
                    for i in range(self.ndims)]

            return ComputeMetaKernel(muls)

        kernels['gradcoru_fpts'] = gradcoru_fpts

        if 'flux' in self.antialias:
            def gradcoru_qpts():
                nupts, nqpts = self.nupts, self.nqpts
                vupts, vqpts = self._vect_upts, self._vect_qpts

                # Exploit the block-diagonal form of the operator
                muls = [self._be.kernel('mul', self.opmat('M7'),
                                        vupts.slice(i*nupts, (i + 1)*nupts),
                                        vqpts.slice(i*nqpts, (i + 1)*nqpts))
                        for i in range(self.ndims)]

                return ComputeMetaKernel(muls)

            kernels['gradcoru_qpts'] = gradcoru_qpts

        # Shock capturing
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')
        if shock_capturing == 'artificial-viscosity':
            tags = {'align'}

            # Register the kernels
            self._be.pointwise.register('pyfr.solvers.baseadvecdiff.kernels.shocksensor')

            # Obtain the scalar variable to be used for shock sensing
            shockvar = 0

            ename = self.basis.name
            weights = get_quadrule(ename, self.cfg.get(f'solver-elements-{ename}', 'soln-pts'), self.nupts).wts
            weights /= np.sum(weights)

            dt_rev = float(self.cfg.get('solver-artificial-viscosity', 'dt_rev', self.cfg.get('solver-time-integrator', 'dt')))
            c_mu = float(self.cfg.get('solver-artificial-viscosity', 'c_mu'))
            cutoff = float(self.cfg.get('solver-artificial-viscosity', 'cutoff', 0.0))
            exp_fac = float(self.cfg.get('solver-artificial-viscosity', 'exp_fac', 1.0))
            vis_coeffs = self.cfg.getliteral('solver-artificial-viscosity', 'vis_coeffs', [1.0]*self.nvars)
            vis_method = self.cfg.get('solver-artificial-viscosity', 'vis_method')

            sensor_type = self.cfg.get('solver-artificial-viscosity', 'sensor', 'rev')
            ubdegs = [sum(dd) for dd in self.basis.ubasis.degrees]

            # Template arguments
            tplargs = dict(
                nvars=self.nvars, nupts=self.nupts, ndims=self.ndims,
                c=self.cfg.items_as('solver-artificial-viscosity', float),
                weights=weights, dt_rev=dt_rev, vis_coeffs=vis_coeffs,
                c_mu=c_mu, order=self.basis.order, vis_method=vis_method,
                cutoff=cutoff, exp_fac=exp_fac, 
                ubdegs=ubdegs, invvdm=self.basis.ubasis.invvdm.T,
                sensor_type=sensor_type
            )

            # Allocate space for the artificial viscosity vector
            self.artvisc = self._be.matrix((self.nupts, self.nvars, self.neles),
                                           extent=nonce + 'artvisc', tags=tags)

            # Apply the sensor to estimate the required artificial viscosity
            kernels['shocksensor'] = lambda: self._be.kernel(
                'shocksensor', tplargs=tplargs, dims=[self.neles],
                du=self._scal_upts_cpy, artvisc=self.artvisc,
                rcpdjac=self.rcpdjac_at('upts')
            )
        else:
            self.artvisc = None

    def get_artvisc_fpts_for_inter(self, eidx, fidx):
        nfp = self.nfacefpts[fidx]
        rmap = self._srtd_face_fpts[fidx][eidx]
        rmap = (0,)*nfp
        cmap = (eidx,)*nfp
        return (self.artvisc.mid,)*nfp, rmap, cmap