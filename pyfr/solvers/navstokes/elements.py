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
            ubdegs = [sum(dd) for dd in self.basis.ubasis.degrees]

            # Get time step
            dt = self.cfg.get('solver-time-integrator', 'dt')

            tensorprodelem = self.basis.nupts == (self.basis.order+1)**self.ndims
            # # If tensor product element, use Legendre modes
            if tensorprodelem:
                [vdm, invvdm] = self.makeLegendreMats()
            # Else use Jacobi modes
            else:
                [vdm, invvdm] = [self.basis.ubasis.vdm.T, self.basis.ubasis.invvdm.T]  

            # DEBUG CODE
            # [vdml, invvdml] = self.makeLegendreMats()  
            # [vdmj, invvdmj] = [self.basis.ubasis.vdm.T, self.basis.ubasis.invvdm.T]      
            # x = self.basis.upts[:,0]  
            # y = self.basis.upts[:,1]
            # f = x
            # q = np.matmul(invvdml, f)
            # q[np.abs(q) < 1e-10] = 0.0
            # q = np.reshape(q, (self.basis.order+1, self.basis.order+1))
            # print(q)
            # q = np.matmul(invvdmj, f)
            # q[np.abs(q) < 1e-10] = 0.0
            # q = np.reshape(q, (self.basis.order+1, self.basis.order+1))
            # print(q)
            # input()
            
            sftplargs = dict(
                nvars=self.nvars, ndims=self.ndims, nupts=self.nupts, svar=shockvar,
                c=self.cfg.items_as('solver-modal-filtering', float),
                order=self.basis.order, ubdegs=ubdegs, srcex=self._src_exprs,
                vdm=vdm, invvdm=invvdm, dt=dt, kexp=self.cfg.get('solver-modal-filtering', 'exponent', '-4'),
                filtermethod=self.cfg.get('solver-modal-filtering', 'filter-method', 'pointwise')
            )


            self.shockcell = self._be.matrix((1, self.neles), tags={'align'}, initval=np.zeros((1, self.neles)))

            #Apply the sensor to estimate the required artificial viscosity
            self.kernels['shocksensorfilter'] = lambda: self._be.kernel(
               'shocksensorfilter', tplargs=sftplargs, dims=[self.neles],
               u=self.scal_upts_inb, shockcell=self.shockcell
            )

            if 'flux' in self.antialias:
                raise ValueError('Cannot use anti-aliasing with modal filtering shock capturing.')

            plocupts = self.ploc_at('upts')

            self.kernels['negdivconf'] = lambda: self._be.kernel(
               'negdivconffilter', tplargs=sftplargs,
               dims=[self.neles], tdivtconf=self.scal_upts_outb,
               rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=self.scal_upts_inb, 
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


    def makeLegendreMats(self):
        upts =  self.basis.upts
        if self.ndims == 2:
            vdm  = np.polynomial.legendre.legvander2d(upts[:,0], upts[:,1], [self.basis.order]*2)
        elif self.ndims == 3:
            vdm  = np.polynomial.legendre.legvander3d(upts[:,0], upts[:,1], upts[:,2], [self.basis.order]*3)
        invvdm = np.linalg.inv(vdm)
        invvdm[np.abs(invvdm) < 1e-7] = 0.0
        return [vdm, invvdm]