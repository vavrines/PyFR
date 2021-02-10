# -*- coding: utf-8 -*-

import numpy as np

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from scipy.optimize import minimize_scalar


class NavierStokesIntInters(BaseAdvectionDiffusionIntInters):
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        # Pointwise template arguments
        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = self.cfg.get('solver', 'viscosity-correction')
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       visc_corr=visc_corr, shock_capturing=shock_capturing,
                       c=self._tpl_c)


        be.pointwise.register('pyfr.solvers.navstokes.kernels.intconu')
        be.pointwise.register('pyfr.solvers.navstokes.kernels.intcflux')

        self.walldist = walldist_at_ploc(self, self._ploc, 'intinters')

        if abs(self._tpl_c['ldg-beta']) == 0.5:
            self.kernels['copy_fpts'] = lambda: ComputeMetaKernel(
                [ele.kernels['_copy_fpts']() for ele in elemap.values()]
            )

        self.kernels['con_u'] = lambda: be.kernel(
            'intconu', tplargs=tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs,
            ulout=self._vect_lhs, urout=self._vect_rhs
        )
        self.kernels['comm_flux'] = lambda: be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            walldist=self.walldist
        )


class NavierStokesMPIInters(BaseAdvectionDiffusionMPIInters):
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        # Pointwise template arguments
        rsolver = cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = cfg.get('solver', 'viscosity-correction')
        shock_capturing = cfg.get('solver', 'shock-capturing')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       visc_corr=visc_corr, shock_capturing=shock_capturing,
                       c=self._tpl_c)

        be.pointwise.register('pyfr.solvers.navstokes.kernels.mpiconu')
        be.pointwise.register('pyfr.solvers.navstokes.kernels.mpicflux')

        self.walldist = walldist_at_ploc(self, self._ploc, 'intinters')

        self.kernels['con_u'] = lambda: be.kernel(
            'mpiconu', tplargs=tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs, ulout=self._vect_lhs
        )
        self.kernels['comm_flux'] = lambda: be.kernel(
            'mpicflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            walldist=self.walldist
        )


class NavierStokesBaseBCInters(BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        # Pointwise template arguments
        rsolver = cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = cfg.get('solver', 'viscosity-correction', 'none')
        shock_capturing = cfg.get('solver', 'shock-capturing', 'none')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                       visc_corr=visc_corr, shock_capturing=shock_capturing,
                       c=self._tpl_c, bctype=self.type,
                       bccfluxstate=self.cflux_state)

        be.pointwise.register('pyfr.solvers.navstokes.kernels.bcconu')
        be.pointwise.register('pyfr.solvers.navstokes.kernels.bccflux')

        self.walldist = walldist_at_ploc(self, self._ploc, 'intinters')

        self.kernels['con_u'] = lambda: be.kernel(
            'bcconu', tplargs=tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, ulout=self._vect_lhs,
            nlin=self._norm_pnorm_lhs, ploc=self._ploc
        )
        self.kernels['comm_flux'] = lambda: be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, gradul=self._vect_lhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            ploc=self._ploc, artviscl=self._artvisc_lhs,
            walldist=self.walldist
        )


class NavierStokesNoSlpIsotWallBCInters(NavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self._tpl_c['cpTw'], = self._eval_opts(['cpTw'])   
        self._tpl_c['wuWall'], = self._eval_opts(['wuWall'])  
        self._tpl_c.update(
            self._exp_opts('uvw'[:self.ndims], lhs,
                           default={'u': 0, 'v': 0, 'w': 0})
        )


class NavierStokesNoSlpAdiaWallBCInters(NavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost'
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self._tpl_c['wuWall'], = self._eval_opts(['wuWall'])  


class NavierStokesSlpAdiaWallBCInters(NavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class NavierStokesCharRiemInvBCInters(NavierStokesBaseBCInters):
    type = 'char-riem-inv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        tplc = self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2] + ['ku','wu'], lhs
        )
        self._tpl_c.update(tplc)


class NavierStokesSupInflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        tplc = self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2] + ['ku','wu'], lhs
        )
        self._tpl_c.update(tplc)


class NavierStokesSupOutflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class NavierStokesSubInflowFrvBCInters(NavierStokesBaseBCInters):
    type = 'sub-in-frv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        tplc = self._exp_opts(
            ['rho', 'u', 'v', 'w'][:self.ndims + 1] + ['ku','wu'], lhs,
            default={'u': 0, 'v': 0, 'w': 0}
        )
        self._tpl_c.update(tplc)


class NavierStokesSubInflowFtpttangBCInters(NavierStokesBaseBCInters):
    type = 'sub-in-ftpttang'
    cflux_state = 'ghost'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        gamma = self.cfg.getfloat('constants', 'gamma')

        # Pass boundary constants to the backend
        self._tpl_c['cpTt'], = self._eval_opts(['cpTt'])
        self._tpl_c['pt'], = self._eval_opts(['pt'])
        self._tpl_c['Rdcp'] = (gamma - 1.0)/gamma

        # Calculate u, v velocity components from the inflow angle
        theta = self._eval_opts(['theta'])[0]*np.pi/180.0
        velcomps = np.array([np.cos(theta), np.sin(theta), 1.0])

        # Adjust u, v and calculate w velocity components for 3-D
        if self.ndims == 3:
            phi = self._eval_opts(['phi'])[0]*np.pi/180.0
            velcomps[:2] *= np.sin(phi)
            velcomps[2] *= np.cos(phi)

        self._tpl_c['vc'] = velcomps[:self.ndims]


class NavierStokesSubOutflowBCInters(NavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self._tpl_c.update(self._exp_opts(['p'], lhs))

def walldist_at_ploc(self, ploc, nonce):
    global x0, y0
    geo = self.cfg.get('solver', 'geometry')
    plocdata = ploc.get()
    walldist = np.zeros_like(plocdata)    

    coords = plocdata.swapaxes(0, 1)
    npts = len(coords[:,0]) # shape of coords x vector to get # pts
    for i in range(npts):
        [x,y,z] = coords[i,:]
        if geo == 'cylinder':
            d = np.sqrt(x**2 + y**2) - 0.5
        elif geo == 'squarecylinder':
            d1 = max(0., abs(x) - 0.5)
            d2 = max(0., abs(y) - 0.5) 
            d = np.sqrt(d1**2 + d2**2)
        elif geo == 'tandsphere':
            d = min(np.sqrt(x**2 + y**2), (x-10.)**2 + y**2) - 0.5
        elif geo == 'cube':
            d1 = max(0., np.abs(x) - 0.5)
            d2 = max(0., np.abs(y) - 0.5)
            d3 = max(0., np.abs(z) - 0.5)
            d = np.sqrt(d1**2 + d2**2 + d3**2)
            d = min(d, y)
        elif geo == 'TGV':
            d = 100000000
        elif geo == 'channel':
            d = abs(y) - 1.0
        elif geo == 'periodichill':
            x0 = x
            y0 = y
            if x >= 54.0/28.0 and x <= 54.0/28.0 + 9.0:
                d = y
            else:
                f = minimize_scalar(periodic_hill_obj, tol=1e-4, method='brent', options={'xtol': 1e-04, 'maxiter': 10})
                d = periodic_hill_obj(f.x)
            d = min(d, 3.035 - y)
            d = max(0.0, d)

        walldist[:,i] = d

    walldist  = self._be.matrix(np.shape(plocdata), tags={'align'}, extent= 'walldist' + nonce, initval=walldist)
    return walldist


def periodic_hill_geo(xx):
    x = xx*28.0
    if x <= 9.:
        return min(28., 2.800000000000E+01 + 0.000000000000E+00*x + 6.775070969851E-03*x**2 -2.124527775800E-03*x**3)
    if x > 9. and x <= 14.:
        return 2.507355893131E+01 +9.754803562315E-01*x -1.016116352781E-01*x**2  +1.889794677828E-03*x**3
    if x > 14. and x <= 20.:
        return 2.579601052357E+01      +8.206693007457E-01*x      -9.055370274339E-02*x**2  +1.626510569859E-03*x**3 
    if x > 20. and x <= 30.:
        return 4.046435022819E+01      -1.379581654948E+00*x      +1.945884504128E-02*x**2  -2.070318932190E-04*x**3
    if x > 30. and x <= 40.:
        return 1.792461334664E+01      +8.743920332081E-01*x      -5.567361123058E-02*x**2  +6.277731764683E-04*x**3
    if x > 40. and x <= 54.:
        return max(0., 5.639011190988E+01      -2.010520359035E+00*x       +1.644919857549E-02*x**2  +2.674976141766E-05*x**3)
    if x > 54. and x <= 306.:
        return 0.0
    if x > 306:
        return periodic_hill_geo(306 +54 - x)

def periodic_hill_obj(x):
    global x0, y0
    y = periodic_hill_geo(x)
    return np.sqrt((x - x0)**2 + (y - y0)**2)
