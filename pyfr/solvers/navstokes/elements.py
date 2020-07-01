# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements
import numpy as np

class NavierStokesElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.tflux')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float))

		# Shock capturing
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')
        if shock_capturing == 'modal-filtering':

            # Register the kernels
            self._be.pointwise.register('pyfr.solvers.navstokes.kernels.shocksensorfilter')
            self._be.pointwise.register('pyfr.solvers.navstokes.kernels.negdivconffilter')

            # Use density for sensing
            shockvar = 0

            # Obtain the degrees of the polynomial modes in the basis
            ubdegs = np.array([sum(dd) for dd in self.basis.ubasis.degrees])

            # Get time step
            dt = self.cfg.get('solver-time-integrator', 'dt')

            # Get Legendre Vandermonde matrix 
            vdm = self.basis.ubasis.vdm.T

            # vdm is currently using orthonormal Legendre modes (int(u_i*u_i) = 1) which have a factor of sqrt(i + 0.5)
            # Transform to orthogonal Legendre modes
            normfac = np.zeros_like(ubdegs, dtype=float)
            for i,dd in enumerate(self.basis.ubasis.degrees):
            	normfac[i] = np.prod(np.sqrt(1./(np.array(dd)+0.5)))

            # Compute orthogonal Vandermonde matrix and inverse
            vdm = normfac*vdm
            invvdm = np.linalg.inv(vdm)


            # Get target decay rates
            decay = self.cfg.get('solver-modal-filtering', 'decay', 'exponential')
            # Exponential rates assume c_inf continuity with a decay rate of 2.0 (since modes are squared)
            if decay == 'exponential':
            	decayfactor = self.cfg.getfloat('solver-modal-filtering', 'decay-factor', 2.0)
            	kexp = self.cfg.getfloat('solver-modal-filtering', 'filter-exponent', 2.0)
            	targetratios = np.exp(-decayfactor*(ubdegs[1:]))
            # Power law rates assume c_0 continuity with a decay rate of 4.0 (since modes are squared)
            elif decay == 'power':
            	decayfactor = self.cfg.getfloat('solver-modal-filtering', 'decay-factor', 4.0)
            	targetratios = (ubdegs[1:]+1.)**(-decayfactor)
            	kexp = None

            # Template arguments
            sftplargs = dict(
                nvars=self.nvars, ndims=self.ndims, nupts=self.nupts, svar=shockvar,
                c=self.cfg.items_as('solver-modal-filtering', float),
                order=self.basis.order, ubdegs=ubdegs, srcex=self._src_exprs,
                vdm=vdm, invvdm=invvdm, dt=dt, targetratios=targetratios, kexp=kexp,
                filtermethod=self.cfg.get('solver-modal-filtering', 'filter-method', 'exponential')
            )

            # Shockcell is an elementwise value (0,1) dictating if there is a shock in the cell (1)
            self.shockcell = self._be.matrix((1, self.neles), tags={'align'}, initval=np.zeros((1, self.neles)))

            # Currently the sensor is always on and applies filter everywhere
            self.kernels['shocksensorfilter'] = lambda: self._be.kernel(
               'shocksensorfilter', tplargs=sftplargs, dims=[self.neles],
               u=self.scal_upts_inb, shockcell=self.shockcell
            )

            # Flux points need to be aligned with solution points
            if 'flux' in self.antialias:
                raise ValueError('Cannot use anti-aliasing with modal filtering shock capturing.')

            # Change negdivconf kernel to include filtering 
            self.kernels['negdivconf'] = lambda: self._be.kernel(
               'negdivconffilter', tplargs=sftplargs,
               dims=[self.neles], tdivtconf=self.scal_upts_outb,
               rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb, 
               shockcell=self.shockcell
            )

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
        else:
            self.kernels['tdisf'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )
