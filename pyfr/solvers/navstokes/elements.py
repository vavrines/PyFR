# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.euler.elements import BaseFluidElements
import numpy as np
from scipy import interpolate

class NavierStokesElements(BaseFluidElements, BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.zero_interior_pressure')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.negdivconf_ns')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.compute_divergence')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.correct_pressure')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        dt = self.cfg.get('solver-time-integrator', 'dt')

        ILM = InvLapMat(self)
        FLM = InterLapMat(self)
        [Mx, My, Mz] = GradMats(self)
        [Gx, Gy, Gz] = InterGradMats(self)
        [Lx, Ly, Lz] = LOGradMats(self)

        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nupts=self.nupts,
                       nfpts=self.nfpts, shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float), srcex=self._src_exprs,
                       dt=dt, ILM=ILM, FLM=FLM, 
                       Mx=Mx, My=My, Mz=Mz, 
                       Gx=Gx, Gy=Gy, Gz=Gz, 
                       Lx=Lx, Ly=Ly, Lz=Lz)

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


        self.kernels['zero_interior_pressure'] = lambda: self._be.kernel(
            'zero_interior_pressure', tplargs=tplargs, dims=[self.nupts, self.neles],
            u=self.scal_upts_inb
        )

        plocupts = self.ploc_at('upts') if self._ploc_in_src_exprs else None
        self.kernels['negdivconf_ns'] = lambda: self._be.kernel(
            'negdivconf_ns', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=self.scal_upts_inb
        )

        self.kernels['compute_divergence'] = lambda: self._be.kernel(
            'compute_divergence', tplargs=tplargs,
            dims=[self.nupts, self.neles], u=self.scal_upts_outb,
            gradu=self._vect_upts
        )

        self.kernels['correct_pressure'] = lambda: self._be.kernel(
            'correct_pressure', tplargs=tplargs,
            dims=[self.neles], uoutb=self.scal_upts_outb,
            uinb=self.scal_upts_inb, rcpdjac=self.rcpdjac_at('upts'),
            ufpts=self._scal_fpts_cpy, ucpy=self._scal_upts_cpy
        )

# Inverse of interior Laplacian matrix
def InvLapMat(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    M = np.zeros((self.nupts, self.nupts)) 

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    if self.ndims == 2:
        n = 0
        uidx = lambda i, j: i + j*(p+1)

        for j in range(p+1):
            for i in range(p+1):
                valsi = np.zeros(p+3)
                valsi[i+1] = 1.0
                dlagi = interpolate.lagrange(dfrpts, valsi).deriv().deriv()

                valsj = np.zeros(p+3)
                valsj[j+1] = 1.0
                dlagj = interpolate.lagrange(dfrpts, valsj).deriv().deriv()

                for ii in range(p+1):
                    M[uidx(ii,j), n] += dlagi(upts_1d[ii])

                for jj in range(p+1):
                    M[uidx(i,jj), n] += dlagj(upts_1d[jj])

                n += 1


    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2


        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+3)
                    valsi[i+1] = 1.0
                    dlagi = interpolate.lagrange(dfrpts, valsi).deriv().deriv()

                    valsj = np.zeros(p+3)
                    valsj[j+1] = 1.0
                    dlagj = interpolate.lagrange(dfrpts, valsj).deriv().deriv()

                    valsk = np.zeros(p+3)
                    valsk[k+1] = 1.0
                    dlagk = interpolate.lagrange(dfrpts, valsk).deriv().deriv()

                    for ii in range(p+1):
                        M[uidx(ii,j,k), n] += dlagi(upts_1d[ii])

                    for jj in range(p+1):
                        M[uidx(i,jj,k), n] += dlagj(upts_1d[jj])

                    for kk in range(p+1):
                        M[uidx(i,j,kk), n] += dlagk(upts_1d[kk])

                    n += 1

    M = np.linalg.inv(M)
    return M

# Interface Laplacian matrix
def InterLapMat(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    M = np.zeros((self.nupts, self.nfpts)) 

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    vals = np.zeros(p+3)
    vals[0] = 1.
    g_left = interpolate.lagrange(dfrpts, vals).deriv().deriv()

    vals = np.zeros(p+3)
    vals[-1] = 1.
    g_right = interpolate.lagrange(dfrpts, vals).deriv().deriv()

    if self.ndims == 2:
        # QUAD FACE ORDERING:
        # 0 -> y = -1
        # 1 -> x =  1
        # 2 -> y =  1
        # 3 -> x = -1
        faceidx = lambda face, idx: face*(p+1) + idx
        uidx = lambda xidx, yidx: xidx + yidx*(p+1)
        for yidx in range(p+1):
            for xidx in range(p+1):
                solidx = uidx(xidx, yidx)
                fxm_idx = faceidx(3,yidx)
                fxp_idx = faceidx(1,yidx)
                fym_idx = faceidx(0,xidx)
                fyp_idx = faceidx(2,xidx)
                M[solidx, fxm_idx] = g_left(upts_1d[xidx])
                M[solidx, fym_idx] = g_left(upts_1d[yidx])
                M[solidx, fxp_idx] = g_right(upts_1d[xidx])
                M[solidx, fyp_idx] = g_right(upts_1d[yidx])

    elif self.ndims == 3:
        # HEX FACE ORDERING:
        # 0 -> z = -1
        # 1 -> y = -1
        # 2 -> x =  1
        # 3 -> y =  1
        # 4 -> x = -1
        # 5 -> z =  1
        faceidx = lambda face, a_idx, b_idx: face*(p+1)**2 + a_idx + b_idx*(p+1)
        uidx = lambda xidx, yidx, zidx: xidx + yidx*(p+1) + zidx*((p+1)**2)
        for zidx in range(p+1):
            for yidx in range(p+1):
                for xidx in range(p+1):
                    solidx = uidx(xidx, yidx, zidx)
                    fxm_idx = faceidx(4,yidx,zidx)
                    fxp_idx = faceidx(2,yidx,zidx)
                    fym_idx = faceidx(1,xidx,zidx)
                    fyp_idx = faceidx(3,xidx,zidx)
                    fzm_idx = faceidx(0,xidx,yidx)
                    fzp_idx = faceidx(5,xidx,yidx)
                    M[solidx, fxm_idx] = g_left(upts_1d[xidx])
                    M[solidx, fym_idx] = g_left(upts_1d[yidx])
                    M[solidx, fzm_idx] = g_left(upts_1d[zidx])
                    M[solidx, fxp_idx] = g_right(upts_1d[xidx])
                    M[solidx, fyp_idx] = g_right(upts_1d[yidx])
                    M[solidx, fzp_idx] = g_right(upts_1d[zidx])
    return M

# Interface gradient matrices
def InterGradMats(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mx = np.zeros((self.nupts, self.nfpts)) 
    My = np.zeros((self.nupts, self.nfpts)) 
    Mz = np.zeros((self.nupts, self.nfpts))

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    vals = np.zeros(p+3)
    vals[0] = 1.
    g_left = interpolate.lagrange(dfrpts, vals).deriv()

    vals = np.zeros(p+3)
    vals[-1] = 1.
    g_right = interpolate.lagrange(dfrpts, vals).deriv()

    if self.ndims == 2:
        # QUAD FACE ORDERING:
        # 0 -> y = -1
        # 1 -> x =  1
        # 2 -> y =  1
        # 3 -> x = -1
        faceidx = lambda face, idx: face*(p+1) + idx
        uidx = lambda xidx, yidx: xidx + yidx*(p+1)
        for yidx in range(p+1):
            for xidx in range(p+1):
                solidx = uidx(xidx, yidx)
                fxm_idx = faceidx(3,yidx)
                fxp_idx = faceidx(1,yidx)
                fym_idx = faceidx(0,xidx)
                fyp_idx = faceidx(2,xidx)
                Mx[solidx, fxm_idx] = g_left(upts_1d[xidx])
                My[solidx, fym_idx] = g_left(upts_1d[yidx])
                Mx[solidx, fxp_idx] = g_right(upts_1d[xidx])
                My[solidx, fyp_idx] = g_right(upts_1d[yidx])

    elif self.ndims == 3:
        # HEX FACE ORDERING:
        # 0 -> z = -1
        # 1 -> y = -1
        # 2 -> x =  1
        # 3 -> y =  1
        # 4 -> x = -1
        # 5 -> z =  1
        faceidx = lambda face, a_idx, b_idx: face*(p+1)**2 + a_idx + b_idx*(p+1)
        uidx = lambda xidx, yidx, zidx: xidx + yidx*(p+1) + zidx*((p+1)**2)
        for zidx in range(p+1):
            for yidx in range(p+1):
                for xidx in range(p+1):
                    solidx = uidx(xidx, yidx, zidx)
                    fxm_idx = faceidx(4,yidx,zidx)
                    fxp_idx = faceidx(2,yidx,zidx)
                    fym_idx = faceidx(1,xidx,zidx)
                    fyp_idx = faceidx(3,xidx,zidx)
                    fzm_idx = faceidx(0,xidx,yidx)
                    fzp_idx = faceidx(5,xidx,yidx)
                    Mx[solidx, fxm_idx] = g_left(upts_1d[xidx])
                    My[solidx, fym_idx] = g_left(upts_1d[yidx])
                    Mz[solidx, fzm_idx] = g_left(upts_1d[zidx])
                    Mx[solidx, fxp_idx] = g_right(upts_1d[xidx])
                    My[solidx, fyp_idx] = g_right(upts_1d[yidx])
                    Mz[solidx, fzp_idx] = g_right(upts_1d[zidx])
    return [Mx, My, Mz]

def LOGradMats(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mx = np.zeros((self.nupts, self.nupts)) 
    My = np.zeros((self.nupts, self.nupts)) 
    Mz = np.zeros((self.nupts, self.nupts)) if self.ndims == 3 else None


    if self.ndims == 2:
        n = 0
        uidx = lambda i, j: i + j*(p+1)

        for j in range(p+1):
            for i in range(p+1):
                valsi = np.zeros(p+1)
                valsi[i] = 1.0
                dlagi = interpolate.lagrange(upts_1d, valsi).deriv()

                valsj = np.zeros(p+1)
                valsj[j] = 1.0
                dlagj = interpolate.lagrange(upts_1d, valsj).deriv()

                for ii in range(p+1):
                    Mx[uidx(ii,j), n] += dlagi(upts_1d[ii])

                for jj in range(p+1):
                    My[uidx(i,jj), n] += dlagj(upts_1d[jj])

                n += 1


    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2


        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+1)
                    valsi[i] = 1.0
                    dlagi = interpolate.lagrange(upts_1d, valsi).deriv()

                    valsj = np.zeros(p+1)
                    valsj[j] = 1.0
                    dlagj = interpolate.lagrange(upts_1d, valsj).deriv()

                    valsk = np.zeros(p+1)
                    valsk[k] = 1.0
                    dlagk = interpolate.lagrange(upts_1d, valsk).deriv()

                    for ii in range(p+1):
                        Mx[uidx(ii,j,k), n] += dlagi(upts_1d[ii])

                    for jj in range(p+1):
                        My[uidx(i,jj,k), n] += dlagj(upts_1d[jj])

                    for kk in range(p+1):
                        Mz[uidx(i,j,kk), n] += dlagk(upts_1d[kk])

                    n += 1

    
    return [Mx, My, Mz]

def GradMats(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mx = np.zeros((self.nupts, self.nupts)) 
    My = np.zeros((self.nupts, self.nupts)) 
    Mz = np.zeros((self.nupts, self.nupts)) if self.ndims == 3 else None

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    if self.ndims == 2:
        n = 0
        uidx = lambda i, j: i + j*(p+1)

        for j in range(p+1):
            for i in range(p+1):
                valsi = np.zeros(p+3)
                valsi[i+1] = 1.0
                dlagi = interpolate.lagrange(dfrpts, valsi).deriv()

                valsj = np.zeros(p+3)
                valsj[j+1] = 1.0
                dlagj = interpolate.lagrange(dfrpts, valsj).deriv()

                for ii in range(p+1):
                    Mx[uidx(ii,j), n] += dlagi(upts_1d[ii])

                for jj in range(p+1):
                    My[uidx(i,jj), n] += dlagj(upts_1d[jj])

                n += 1


    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2


        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+3)
                    valsi[i+1] = 1.0
                    dlagi = interpolate.lagrange(dfrpts, valsi).deriv()

                    valsj = np.zeros(p+3)
                    valsj[j+1] = 1.0
                    dlagj = interpolate.lagrange(dfrpts, valsj).deriv()

                    valsk = np.zeros(p+3)
                    valsk[k+1] = 1.0
                    dlagk = interpolate.lagrange(dfrpts, valsk).deriv()

                    for ii in range(p+1):
                        Mx[uidx(ii,j,k), n] += dlagi(upts_1d[ii])

                    for jj in range(p+1):
                        My[uidx(i,jj,k), n] += dlagj(upts_1d[jj])

                    for kk in range(p+1):
                        Mz[uidx(i,j,kk), n] += dlagk(upts_1d[kk])

                    n += 1

    
    return [Mx, My, Mz]

        