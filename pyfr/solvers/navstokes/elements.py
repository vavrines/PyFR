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
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.negdivconf_inv')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.compute_divergence')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.correct_pressure')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.clean_divergence')

        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        dt = self.cfg.get('solver-time-integrator', 'dt')

        [Dx, Dy, Dz] = DFRGradMat(self)
        [Sx, Sy, Sz] = SolGradMat(self)
        [Fx, Fy, Fz] = FluxGradMat(self)
        V = ProjMat(self)

        self.invLapMat = self._be.const_matrix(self.makeSolLapMats(Dx, Dy, Dz, Sx, Sy, Sz), tags={'align'})
        self.fluxLapMat = self._be.const_matrix(self.makeFluxLapMats(Dx, Dy, Dz, Fx, Fy, Fz), tags={'align'})

        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nupts=self.nupts,
                       nfpts=self.nfpts, shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float), srcex=self._src_exprs,
                       dt=dt, V=V)

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
        self.kernels['negdivconf_inv'] = lambda: self._be.kernel(
            'negdivconf_inv', tplargs=tplargs,
            dims=[self.nupts, self.neles], uout=self.scal_upts_outb,
            uin=self.scal_upts_inb
        )

        self.kernels['compute_divergence'] = lambda: self._be.kernel(
            'compute_divergence', tplargs=tplargs,
            dims=[self.nupts, self.neles], u=self.scal_upts_outb,
            gradu=self._vect_upts
        )

        self.kernels['correct_pressure'] = lambda: self._be.kernel(
            'correct_pressure', tplargs=tplargs,
            dims=[self.neles], uoutb=self.scal_upts_outb,
            uinb=self.scal_upts_inb, ufpts=self._scal_fpts_cpy, 
            ILM=self.invLapMat, FLM=self.fluxLapMat
        )

        self.kernels['clean_divergence'] = lambda: self._be.kernel(
            'clean_divergence', tplargs=tplargs,
            dims=[self.neles], uoutb=self.scal_upts_outb
        )

    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        return grad_cons

    def makeSolLapMats(self, Dx, Dy, Dz, Sx, Sy, Sz):
        urcpdjac = self.rcpdjac_at_np('upts') # (nupts, nelems)
        frcpdjac = self.rcpdjac_at_np('fpts') # (nfpts, nelems)
        drcpdjac = np.concatenate((urcpdjac, frcpdjac), axis=0)
        usmats = self.smat_at_np('upts') # (ndims, nupts, ndims, nelems)
        fsmats = self.smat_at_np('fpts') # (ndims, nfpts, ndims, nelems)
        dsmats = np.concatenate((usmats, fsmats), axis=1)

        [_, nelems] = np.shape(urcpdjac)
        invLapMat = np.zeros((self.nupts, self.nupts, nelems))


        uniform_grid = self.cfg.get('solver', 'grid-type', None) == 'uniform'

        if self.ndims == 2:
            for eidx in range(nelems):
                if uniform_grid and eidx > 0:
                    invLapMat[:,:,eidx] = invLapMat[:,:,0]
                else:
                    sxx = usmats[0, :, 0, eidx]
                    sxy = usmats[1, :, 0, eidx]
                    syx = usmats[0, :, 1, eidx]
                    syy = usmats[1, :, 1, eidx]

                    Mx = (drcpdjac[:, eidx]*(sxx*Sx + sxy*Sy).T).T
                    My = (drcpdjac[:, eidx]*(syx*Sx + syy*Sy).T).T

                    sxx = dsmats[0, :, 0, eidx]
                    sxy = dsmats[1, :, 0, eidx]
                    syx = dsmats[0, :, 1, eidx]
                    syy = dsmats[1, :, 1, eidx]


                    Mxx = (urcpdjac[:, eidx]*(sxx*Dx + sxy*Dy).T).T @ Mx
                    Myy = (urcpdjac[:, eidx]*(syx*Dx + syy*Dy).T).T @ My

                    M = Mxx + Myy
                    invLapMat[:,:,eidx] = np.linalg.inv(M)
        elif self.ndims == 3:
            for eidx in range(nelems):
                if uniform_grid and eidx > 0:
                    invLapMat[:,:,eidx] = invLapMat[:,:,0]
                else:
                    sxx = usmats[0, :, 0, eidx]
                    sxy = usmats[1, :, 0, eidx]
                    sxz = usmats[2, :, 0, eidx]
                    syx = usmats[0, :, 1, eidx]
                    syy = usmats[1, :, 1, eidx]
                    syz = usmats[2, :, 1, eidx]
                    szx = usmats[0, :, 2, eidx]
                    szy = usmats[1, :, 2, eidx]
                    szz = usmats[2, :, 2, eidx]

                    Mx = (drcpdjac[:, eidx]*(sxx*Sx + sxy*Sy + sxz*Sz).T).T
                    My = (drcpdjac[:, eidx]*(syx*Sx + syy*Sy + syz*Sz).T).T
                    Mz = (drcpdjac[:, eidx]*(szx*Sx + szy*Sy + szz*Sz).T).T

                    sxx = dsmats[0, :, 0, eidx]
                    sxy = dsmats[1, :, 0, eidx]
                    sxz = dsmats[2, :, 0, eidx]
                    syx = dsmats[0, :, 1, eidx]
                    syy = dsmats[1, :, 1, eidx]
                    syz = dsmats[2, :, 1, eidx]
                    szx = dsmats[0, :, 2, eidx]
                    szy = dsmats[1, :, 2, eidx]
                    szz = dsmats[2, :, 2, eidx]


                    Mxx = (urcpdjac[:, eidx]*(sxx*Dx + sxy*Dy + sxz*Dz).T).T @ Mx
                    Myy = (urcpdjac[:, eidx]*(syx*Dx + syy*Dy + syz*Dz).T).T @ My
                    Mzz = (urcpdjac[:, eidx]*(szx*Dx + szy*Dy + szz*Dz).T).T @ Mz

                    M = Mxx + Myy + Mzz
                    invLapMat[:,:,eidx] = np.linalg.inv(M)
        return invLapMat

    def makeFluxLapMats(self, Dx, Dy, Dz, Fx, Fy, Fz):
        urcpdjac = self.rcpdjac_at_np('upts') # (nupts, nelems)
        frcpdjac = self.rcpdjac_at_np('fpts') # (nfpts, nelems)
        drcpdjac = np.concatenate((urcpdjac, frcpdjac), axis=0)
        usmats = self.smat_at_np('upts') # (ndims, nupts, ndims, nelems)
        fsmats = self.smat_at_np('fpts') # (ndims, nfpts, ndims, nelems)
        dsmats = np.concatenate((usmats, fsmats), axis=1)

        [_, nelems] = np.shape(urcpdjac)
        LapMat = np.zeros((self.nupts, self.nfpts, nelems))


        uniform_grid = self.cfg.get('solver', 'grid-type', None) == 'uniform'

        if self.ndims == 2:
            for eidx in range(nelems):
                if uniform_grid and eidx > 0:
                    LapMat[:,:,eidx] = LapMat[:,:,0]
                else:
                    sxx = fsmats[0, :, 0, eidx]
                    sxy = fsmats[1, :, 0, eidx]
                    syx = fsmats[0, :, 1, eidx]
                    syy = fsmats[1, :, 1, eidx]

                    Mx = (drcpdjac[:, eidx]*(sxx*Fx + sxy*Fy).T).T
                    My = (drcpdjac[:, eidx]*(syx*Fx + syy*Fy).T).T

                    sxx = dsmats[0, :, 0, eidx]
                    sxy = dsmats[1, :, 0, eidx]
                    syx = dsmats[0, :, 1, eidx]
                    syy = dsmats[1, :, 1, eidx]


                    Mxx = (urcpdjac[:, eidx]*(sxx*Dx + sxy*Dy).T).T @ Mx
                    Myy = (urcpdjac[:, eidx]*(syx*Dx + syy*Dy).T).T @ My

                    LapMat[:,:,eidx] = Mxx + Myy
        elif self.ndims == 3:
            for eidx in range(nelems):
                if uniform_grid and eidx > 0:
                    LapMat[:,:,eidx] = LapMat[:,:,0]
                else:
                    sxx = fsmats[0, :, 0, eidx]
                    sxy = fsmats[1, :, 0, eidx]
                    sxz = fsmats[2, :, 0, eidx]
                    syx = fsmats[0, :, 1, eidx]
                    syy = fsmats[1, :, 1, eidx]
                    syz = fsmats[2, :, 1, eidx]
                    szx = fsmats[0, :, 2, eidx]
                    szy = fsmats[1, :, 2, eidx]
                    szz = fsmats[2, :, 2, eidx]

                    Mx = (drcpdjac[:, eidx]*(sxx*Fx + sxy*Fy + sxz*Fz).T).T
                    My = (drcpdjac[:, eidx]*(syx*Fx + syy*Fy + syz*Fz).T).T
                    Mz = (drcpdjac[:, eidx]*(szx*Fx + szy*Fy + szz*Fz).T).T

                    sxx = dsmats[0, :, 0, eidx]
                    sxy = dsmats[1, :, 0, eidx]
                    sxz = dsmats[2, :, 0, eidx]
                    syx = dsmats[0, :, 1, eidx]
                    syy = dsmats[1, :, 1, eidx]
                    syz = dsmats[2, :, 1, eidx]
                    szx = dsmats[0, :, 2, eidx]
                    szy = dsmats[1, :, 2, eidx]
                    szz = dsmats[2, :, 2, eidx]


                    Mxx = (urcpdjac[:, eidx]*(sxx*Dx + sxy*Dy + sxz*Dz).T).T @ Mx
                    Myy = (urcpdjac[:, eidx]*(syx*Dx + syy*Dy + syz*Dz).T).T @ My
                    Mzz = (urcpdjac[:, eidx]*(szx*Dx + szy*Dy + szz*Dz).T).T @ Mz

                    LapMat[:,:,eidx] = Mxx + Myy + Mzz
        return LapMat


def ProjMat(self):
    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    upts = self.basis.ubasis.pts

    def mat(x,y,z):
        if p == 2:
            out = [ [1, 0],
                    [3*y, 0],     
                    [0.5*(15*y**2) - 1.5, 0],     
                    [0, -1],     
                    [x, -y],     
                    [3*x*y, 0.5 - 0.5*(3*y**2)],    
                    [0, -3*x],     
                    [0.5*(3*x**2) - 0.5, -3*x*y],     
                    [0, 1.5 - 0.5*(15*x**2)] ]
        elif p == 3:
            out = [ [                     1,                      0],
                    [                   3*y,                      0],
                    [      (15*y**2)/2.0 - 3/2.0,                      0],
                    [ (35*y**3)/2.0 - (15*y)/2.0,                      0],
                    [                     0,                     -1],
                    [                     x,                     -y],
                    [                 3*x*y,        1/2.0 - (3*y**2)/2.0],
                    [  x*((15*y**2)/2.0 - 3/2.0),    (3*y)/2.0 - (5*y**3)/2.0],
                    [                     0,                   -3*x],
                    [       (3*x**2)/2.0 - 1/2.0,                 -3*x*y],
                    [ 3*y*((3*x**2)/2.0 - 1/2.0), -3*x*((3*y**2)/2.0 - 1/2.0)],
                    [                     0,       3/2.0 - (15*x**2)/2.0],
                    [   (5*x**3)/2.0 - (3*x)/2.0,  -y*((15*x**2)/2.0 - 3/2.0)],
                    [                     0,  (15*x)/2.0 - (35*x**3)/2.0]]
        elif p == 4:
            out = [ [                                    1,                                     0],
                    [                                  3*y,                                     0],
                    [                     (15*y**2)/2.0 - 3/2.0,                                     0],
                    [                (35*y**3)/2.0 - (15*y)/2.0,                                     0],
                    [     (315*y**4)/8.0 - (105*y**2)/4.0 + 15/8.0,                                     0],
                    [                                    0,                                    -1],
                    [                                    x,                                    -y],
                    [                                3*x*y,                       1/2.0 - (3*y**2)/2.0],
                    [                 x*((15*y**2)/2.0 - 3/2.0),                   (3*y)/2.0 - (5*y**3)/2.0],
                    [           -x*((15*y)/2.0 - (35*y**3)/2.0),       - (35*y**4)/8.0 + (15*y**2)/4.0 - 3/8.0],
                    [                                    0,                                  -3*x],
                    [                      (3*x**2)/2.0 - 1/2.0,                                -3*x*y],
                    [                3*y*((3*x**2)/2.0 - 1/2.0),                -3*x*((3*y**2)/2.0 - 1/2.0)],
                    [ ((3*x**2)/2.0 - 1/2.0)*((15*y**2)/2.0 - 3/2.0),             3*x*((3*y)/2.0 - (5*y**3)/2.0)],
                    [                                    0,                      3/2.0 - (15*x**2)/2.0],
                    [                  (5*x**3)/2.0 - (3*x)/2.0,                 -y*((15*x**2)/2.0 - 3/2.0)],
                    [           -3*y*((3*x)/2.0 - (5*x**3)/2.0), -((15*x**2)/2.0 - 3/2.0)*((3*y**2)/2.0 - 1/2.0)],
                    [                                    0,                 (15*x)/2.0 - (35*x**3)/2.0],
                    [        (35*x**4)/8.0 - (15*x**2)/4.0 + 3/8.0,             y*((15*x)/2.0 - (35*x**3)/2.0)],
                    [                                    0,    - (315*x**4)/8.0 + (105*x**2)/4.0 - 15/8.0]]
        elif p == 5:
            out = [ [                                          1,                                         0],
                    [                                        3*y,                                         0],
                    [                           (15*y**2)/2.0 - 3/2.0,                                         0],
                    [                      (35*y**3)/2.0 - (15*y)/2.0,                                         0],
                    [           (315*y**4)/8.0 - (105*y**2)/4.0 + 15/8.0,                                         0],
                    [      (693*y**5)/8.0 - (315*y**3)/4.0 + (105*y)/8.0,                                         0],
                    [                                          0,                                        -1],
                    [                                          x,                                        -y],
                    [                                      3*x*y,                           1/2.0 - (3*y**2)/2.0],
                    [                       x*((15*y**2)/2.0 - 3/2.0),                       (3*y)/2.0 - (5*y**3)/2.0],
                    [                 -x*((15*y)/2.0 - (35*y**3)/2.0),           - (35*y**4)/8.0 + (15*y**2)/4.0 - 3/8.0],
                    [       x*((315*y**4)/8.0 - (105*y**2)/4.0 + 15/8.0),      - (63*y**5)/8.0 + (35*y**3)/4.0 - (15*y)/8.0],
                    [                                          0,                                      -3*x],
                    [                            (3*x**2)/2.0 - 1/2.0,                                    -3*x*y],
                    [                      3*y*((3*x**2)/2.0 - 1/2.0),                    -3*x*((3*y**2)/2.0 - 1/2.0)],
                    [       ((3*x**2)/2.0 - 1/2.0)*((15*y**2)/2.0 - 3/2.0),                 3*x*((3*y)/2.0 - (5*y**3)/2.0)],
                    [ -((15*y)/2.0 - (35*y**3)/2.0)*((3*x**2)/2.0 - 1/2.0),      -3*x*((35*y**4)/8.0 - (15*y**2)/4.0 + 3/8.0)],
                    [                                          0,                          3/2.0 - (15*x**2)/2.0],
                    [                        (5*x**3)/2.0 - (3*x)/2.0,                     -y*((15*x**2)/2.0 - 3/2.0)],
                    [                 -3*y*((3*x)/2.0 - (5*x**3)/2.0),     -((15*x**2)/2.0 - 3/2.0)*((3*y**2)/2.0 - 1/2.0)],
                    [  -((3*x)/2.0 - (5*x**3)/2.0)*((15*y**2)/2.0 - 3/2.0),  ((3*y)/2.0 - (5*y**3)/2.0)*((15*x**2)/2.0 - 3/2.0)],
                    [                                          0,                     (15*x)/2.0 - (35*x**3)/2.0],
                    [              (35*x**4)/8.0 - (15*x**2)/4.0 + 3/8.0,                 y*((15*x)/2.0 - (35*x**3)/2.0)],
                    [        3*y*((35*x**4)/8.0 - (15*x**2)/4.0 + 3/8.0), ((15*x)/2.0 - (35*x**3)/2.0)*((3*y**2)/2.0 - 1/2.0)],
                    [                                          0,        - (315*x**4)/8.0 + (105*x**2)/4.0 - 15/8.0],
                    [         (63*x**5)/8.0 - (35*x**3)/4.0 + (15*x)/8.0,     -y*((315*x**4)/8.0 - (105*x**2)/4.0 + 15/8.0)],
                    [                                          0,   - (693*x**5)/8.0 + (315*x**3)/4.0 - (105*x)/8.0]]
        return np.array(out)
    
    pmat = mat(upts[:,0], upts[:,1], None)
    [ndpts, _] = np.shape(pmat)
    Vx = np.zeros(((p+1)**self.ndims, ndpts))
    Vy = np.zeros(((p+1)**self.ndims, ndpts))
    Vz = np.zeros(((p+1)**self.ndims, ndpts))
    V = np.zeros((self.ndims*(p+1)**self.ndims, ndpts))

    for i in range(ndpts):
        Vx[:, i] = pmat[i,0]
        Vy[:, i] = pmat[i,1]
    
    if self.ndims == 2:
        V[:(p+1)**self.ndims, :] = Vx
        V[(p+1)**self.ndims:, :] = Vy
    else:
        pass


    Vv = V @ np.linalg.inv(V.T @ V) @ V.T

    return Vv



# Compute gradients from points defined on DFR points (upts + fpts)
def DFRGradMat(self):
    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mx = np.zeros((self.nupts, self.nupts+self.nfpts)) 
    My = np.zeros((self.nupts, self.nupts+self.nfpts)) 
    Mz = np.zeros((self.nupts, self.nupts+self.nfpts)) if self.ndims == 3 else None

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    if self.ndims == 2:
        n = 0
        uidx = lambda i, j: i + j*(p+1)

        # Interior solution points
        for j in range(p+1):
            for i in range(p+1):
                valsi = np.zeros(p+3)
                valsi[i+1] = 1.0
                lagi = interpolate.lagrange(dfrpts, valsi)
                dlagi = lagi.deriv()

                valsj = np.zeros(p+3)
                valsj[j+1] = 1.0
                lagj = interpolate.lagrange(dfrpts, valsj)
                dlagj = lagj.deriv()

                for q in range(p+1):
                    for r in range(p+1):
                        Mx[uidx(q,r), n] += dlagi(upts_1d[q])*lagj(upts_1d[r])
                        My[uidx(q,r), n] += lagi(upts_1d[q])*dlagj(upts_1d[r])

                n += 1
        # Flux points

            # QUAD FACE ORDERING:
            # 0 -> y = -1
            # 1 -> x =  1
            # 2 -> y =  1
            # 3 -> x = -1
        for face in range(4):
            for idx in range(p+1):
                # If negative face
                if face in [0,3]:
                    vals = np.zeros(p+3)
                    vals[0] = 1.0
                # If positive face
                elif face in [1,2]:
                    vals = np.zeros(p+3)
                    vals[-1] = 1.0

                lag = interpolate.lagrange(dfrpts, vals)
                dlag = lag.deriv()


                for q in range(p+1):
                    for r in range(p+1):
                        # If x-faces
                        if face in [1,3]:
                            Mx[uidx(q,r), n] += dlag(upts_1d[q]) if r == idx else 0
                        # If y-faces
                        elif face in [0,2]:
                            My[uidx(q,r), n] += dlag(upts_1d[r]) if q == idx else 0

                n += 1

    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2

        # Interior solution points
        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+3)
                    valsi[i+1] = 1.0
                    lagi = interpolate.lagrange(dfrpts, valsi)
                    dlagi = lagi.deriv()

                    valsj = np.zeros(p+3)
                    valsj[j+1] = 1.0
                    lagj = interpolate.lagrange(dfrpts, valsj)
                    dlagj = lagj.deriv()

                    valsk = np.zeros(p+3)
                    valsk[k+1] = 1.0
                    lagk = interpolate.lagrange(dfrpts, valsk)
                    dlagk = lagk.deriv()

                    for q in range(p+1):
                        for r in range(p+1):
                            for s in range(p+1):
                                Mx[uidx(q,r,s), n] += dlagi(upts_1d[q])*lagj(upts_1d[r])*lagk(upts_1d[s])
                                My[uidx(q,r,s), n] += lagi(upts_1d[q])*dlagj(upts_1d[r])*lagk(upts_1d[s])
                                Mz[uidx(q,r,s), n] += lagi(upts_1d[q])*lagj(upts_1d[r])*dlagk(upts_1d[s])

                    n += 1
        # Flux points

            # HEX FACE ORDERING:
            # 0 -> z = -1
            # 1 -> y = -1
            # 2 -> x =  1
            # 3 -> y =  1
            # 4 -> x = -1
            # 5 -> z =  1

        for face in range(6):
            for bidx in range(p+1):
                for aidx in range(p+1):
                    # If negative face
                    if face in [0,1,4]:
                        vals = np.zeros(p+3)
                        vals[0] = 1.0
                    # If positive face
                    elif face in [2,3,5]:
                        vals = np.zeros(p+3)
                        vals[-1] = 1.0

                    lag = interpolate.lagrange(dfrpts, vals)
                    dlag = lag.deriv()

                    for q in range(p+1):
                        for r in range(p+1):
                            for s in range(p+1):
                                # If x-faces
                                if face in [2,4]:
                                    Mx[uidx(q,r,s), n] += dlag(upts_1d[q]) if r == aidx  and s == bidx else 0
                                # If y-faces
                                elif face in [1,3]:
                                    My[uidx(q,r,s), n] += dlag(upts_1d[r]) if q == aidx  and s == bidx else 0
                                # If z-faces
                                elif face in [0,5]:
                                    Mz[uidx(q,r,s), n] += dlag(upts_1d[s]) if q == aidx  and r == bidx else 0

                    n += 1
    return [Mx, My, Mz]


# Compute gradients from points defined on upts (set fpts = 0)
def SolGradMat(self):
    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    fpts = self.basis.fpts
    Mx = np.zeros((self.nupts+self.nfpts, self.nupts)) 
    My = np.zeros((self.nupts+self.nfpts, self.nupts)) 
    Mz = np.zeros((self.nupts+self.nfpts, self.nupts)) if self.ndims == 3 else None

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d


    if self.ndims == 2:
        n = 0
        uidx = lambda i, j: i + j*(p+1)
        faceidx = lambda face, idx: face*(p+1) + idx

        for j in range(p+1):
            for i in range(p+1):
                valsi = np.zeros(p+3)
                valsi[i+1] = 1.0
                lagi = interpolate.lagrange(dfrpts, valsi)
                dlagi = lagi.deriv()

                valsj = np.zeros(p+3)
                valsj[j+1] = 1.0
                lagj = interpolate.lagrange(dfrpts, valsj)
                dlagj = lagj.deriv()

                for q in range(p+1):
                    for r in range(p+1):
                        Mx[uidx(q,r), n] += dlagi(upts_1d[q])*lagj(upts_1d[r])
                        My[uidx(q,r), n] += lagi(upts_1d[q])*dlagj(upts_1d[r])

                for face in range(4):
                    for idx in range(p+1):
                        [x,y] = fpts[faceidx(face,idx), :]

                        Mx[self.nupts + faceidx(face, idx), n] += dlagi(x)*lagj(y)
                        My[self.nupts + faceidx(face, idx), n] += lagi(x)*dlagj(y)

                n += 1


    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2
        faceidx = lambda face, a_idx, b_idx: face*(p+1)**2 + a_idx + b_idx*(p+1)


        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+3)
                    valsi[i+1] = 1.0
                    lagi = interpolate.lagrange(dfrpts, valsi)
                    dlagi = lagi.deriv()

                    valsj = np.zeros(p+3)
                    valsj[j+1] = 1.0
                    lagj = interpolate.lagrange(dfrpts, valsj)
                    dlagj = lagj.deriv()

                    valsk = np.zeros(p+3)
                    valsk[k+1] = 1.0
                    lagk = interpolate.lagrange(dfrpts, valsk)
                    dlagk = lagk.deriv()


                    for q in range(p+1):
                        for r in range(p+1):
                            for s in range(p+1):
                                Mx[uidx(q,r,s), n] += dlagi(upts_1d[q])*lagj(upts_1d[r])*lagk(upts_1d[s])
                                My[uidx(q,r,s), n] += lagi(upts_1d[q])*dlagj(upts_1d[r])*lagk(upts_1d[s])
                                Mz[uidx(q,r,s), n] += lagi(upts_1d[q])*lagj(upts_1d[r])*dlagk(upts_1d[s])

                    for face in range(6):
                        for bidx in range(p+1):
                            for aidx in range(p+1):
                                [x,y,z] = fpts[faceidx(face,aidx,bidx), :]

                                Mx[self.nupts + faceidx(face, aidx, bidx), n] += dlagi(x)*lagj(y)*lagk(z)
                                My[self.nupts + faceidx(face, aidx, bidx), n] += lagi(x)*dlagj(y)*lagk(z)
                                Mz[self.nupts + faceidx(face, aidx, bidx), n] += lagi(x)*lagj(y)*dlagk(z)
                    n += 1

    return [Mx, My, Mz]


# Compute gradients from points defined on fpts (set upts = 0)
def FluxGradMat(self):
    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack

    Mx = np.zeros((self.nupts+self.nfpts, self.nfpts)) 
    My = np.zeros((self.nupts+self.nfpts, self.nfpts)) 
    Mz = np.zeros((self.nupts+self.nfpts, self.nfpts)) if self.ndims == 3 else None

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    if self.ndims == 2:
        n = 0
        faceidx = lambda face, idx: face*(p+1) + idx
        uidx = lambda xidx, yidx: xidx + yidx*(p+1)

        # QUAD FACE ORDERING:
        # 0 -> y = -1
        # 1 -> x =  1
        # 2 -> y =  1
        # 3 -> x = -1
        for face in range(4):
            for idx in range(p+1):
                # If negative face
                if face in [0,3]:
                    valsa = np.zeros(p+3)
                    valsa[0] = 1.0
                # If positive face
                elif face in [1,2]:
                    valsa = np.zeros(p+3)
                    valsa[-1] = 1.0

                laga = interpolate.lagrange(dfrpts, valsa)
                dlaga = laga.deriv()

                valsb = np.zeros(p+1)
                valsb[idx] = 1.0
                lagb = interpolate.lagrange(upts_1d, valsb)
                dlagb = lagb.deriv()

                for q in range(p+1):
                    for r in range(p+1):
                        # If x-faces
                        if face in [1,3]:
                            Mx[uidx(q,r), n] += dlaga(upts_1d[q]) if r == idx else 0
                        # If y-faces
                        elif face in [0,2]:
                            My[uidx(q,r), n] += dlaga(upts_1d[r]) if q == idx else 0
                for qidx in range(p+1):
                    # If x-faces
                    if face in [1,3]:
                        Mx[self.nupts + faceidx(3, qidx), n] += dlaga(-1)*lagb(upts_1d[qidx])
                        Mx[self.nupts + faceidx(1, qidx), n] += dlaga( 1)*lagb(upts_1d[qidx])
                        My[self.nupts + faceidx(face, qidx), n] += dlagb(upts_1d[qidx])
                    # If y-faces
                    elif face in [0,2]:
                        My[self.nupts + faceidx(0, qidx), n] += dlaga(-1)*lagb(upts_1d[qidx])
                        My[self.nupts + faceidx(2, qidx), n] += dlaga( 1)*lagb(upts_1d[qidx])
                        Mx[self.nupts + faceidx(face, qidx), n] += dlagb(upts_1d[qidx])

                n += 1

    elif self.ndims == 3:
        n = 0
        faceidx = lambda face, a_idx, b_idx: face*(p+1)**2 + a_idx + b_idx*(p+1)
        uidx = lambda xidx, yidx, zidx: xidx + yidx*(p+1) + zidx*((p+1)**2)

        # HEX FACE ORDERING:
        # 0 -> z = -1
        # 1 -> y = -1
        # 2 -> x =  1
        # 3 -> y =  1
        # 4 -> x = -1
        # 5 -> z =  1

        for face in range(6):
            for bidx in range(p+1):
                for aidx in range(p+1):
                    # If negative face
                    if face in [0,1,4]:
                        valsa = np.zeros(p+3)
                        valsa[0] = 1.0
                    # If positive face
                    elif face in [2,3,5]:
                        valsa = np.zeros(p+3)
                        valsa[-1] = 1.0

                    laga = interpolate.lagrange(dfrpts, valsa)
                    dlaga = laga.deriv()

                    valsb = np.zeros(p+1)
                    valsb[aidx] = 1.0
                    lagb = interpolate.lagrange(upts_1d, valsb)
                    dlagb = lagb.deriv()

                    valsc = np.zeros(p+1)
                    valsc[bidx] = 1.0
                    lagc = interpolate.lagrange(upts_1d, valsc)
                    dlagc = lagc.deriv()


                    for q in range(p+1):
                        for r in range(p+1):
                            for s in range(p+1):
                                # If x-faces
                                if face in [2,4]:
                                    Mx[uidx(q,r,s), n] += dlaga(upts_1d[q]) if r == aidx  and s == bidx else 0
                                # If y-faces
                                elif face in [1,3]:
                                    My[uidx(q,r,s), n] += dlaga(upts_1d[r]) if q == aidx  and s == bidx else 0
                                # If z-faces
                                elif face in [0,5]:
                                    Mz[uidx(q,r,s), n] += dlaga(upts_1d[s]) if q == aidx  and r == bidx else 0

                    for qidx in range(p+1):
                        for ridx in range(p+1):
                            # If x-faces
                            if face in [2,4]:
                                Mx[self.nupts + faceidx(4, qidx, ridx), n] += dlaga(-1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                Mx[self.nupts + faceidx(2, qidx, ridx), n] += dlaga( 1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                My[self.nupts + faceidx(face, qidx, ridx), n] += dlagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                Mz[self.nupts + faceidx(face, qidx, ridx), n] += lagb(upts_1d[qidx])*dlagc(upts_1d[ridx])
                            # If y-faces
                            elif face in [1,3]:
                                Mx[self.nupts + faceidx(face, qidx, ridx), n] += dlagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                My[self.nupts + faceidx(1, qidx, ridx), n] += dlaga(-1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                My[self.nupts + faceidx(3, qidx, ridx), n] += dlaga( 1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                Mz[self.nupts + faceidx(face, qidx, ridx), n] += lagb(upts_1d[qidx])*dlagc(upts_1d[ridx])
                            # If z-faces
                            elif face in [0,5]:
                                Mx[self.nupts + faceidx(face, qidx, ridx), n] += dlagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                My[self.nupts + faceidx(face, qidx, ridx), n] += lagb(upts_1d[qidx])*dlagc(upts_1d[ridx])
                                Mz[self.nupts + faceidx(0, qidx, ridx), n] += dlaga(-1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])
                                Mz[self.nupts + faceidx(5, qidx, ridx), n] += dlaga( 1)*lagb(upts_1d[qidx])*lagc(upts_1d[ridx])

                    n += 1
                        
    return [Mx, My, Mz]

