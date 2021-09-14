# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements


class MHDElements(BaseFluidElements, BaseAdvectionDiffusionElements):
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

    @property
    def _soln_in_src_exprs(self):
        return True

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tflux_inv')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tflux_vis')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.get_du')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.negdivconf_mhd')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.negdivconf_f')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.negdivconf_b')
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.enforce_positivity')


        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        dt_rev = float(self.cfg.get('solver-artificial-viscosity', 'dt_rev', self.cfg.get('solver-time-integrator', 'dt')))
        mu_max = float(self.cfg.get('solver-artificial-viscosity', 'mu_max'))
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float),
                       srcex=self._src_exprs, dt_rev=dt_rev,
                       mu_max=mu_max)

        if 'flux' in self.antialias:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux_inv', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux_vis', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, artvisc=self.artvisc
            )
        else:
            self.kernels['tdisf_inv'] = lambda: self._be.kernel(
                'tflux_inv', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )
            self.kernels['tdisf_vis'] = lambda: self._be.kernel(
                'tflux_vis', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, artvisc=self.artvisc
            )


        self.kernels['negdivconf'] = lambda: self._be.kernel(
            'negdivconf_mhd', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb
        )

        self.kernels['negdivconf_f'] = lambda: self._be.kernel(
            'negdivconf_f', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb
        )
        self.kernels['negdivconf_b'] = lambda: self._be.kernel(
            'negdivconf_b', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=self.ploc_at('upts'), u=self.scal_upts_inb
        )

        self.kernels['get_du'] = lambda: self._be.kernel(
            'get_du', tplargs=tplargs,
            dims=[self.nupts, self.neles], u_orig=self._scal_upts_cpy, 
            u_rev=self.scal_upts_inb
        )

        self.kernels['enforce_positivity_interface'] = lambda: self._be.kernel(
            'enforce_positivity', tplargs=tplargs,
            dims=[self.nupts, self.neles], u=self._scal_fpts
        )

        self.kernels['enforce_positivity_interior'] = lambda: self._be.kernel(
            'enforce_positivity', tplargs=tplargs,
            dims=[self.nupts, self.neles], u=self.scal_upts_inb
        )