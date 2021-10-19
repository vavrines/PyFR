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

        refsollapmats = refSolLapMat(self)
        reffluxlapmats = refFluxLapMat(self)
        self.invLapMat = self._be.const_matrix(self.makeSolLapMats(refsollapmats), tags={'align'})
        self.fluxLapMat = self._be.const_matrix(self.makeFluxLapMats(reffluxlapmats), tags={'align'})

        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nupts=self.nupts,
                       nfpts=self.nfpts, shock_capturing=shock_capturing, visc_corr=visc_corr,
                       c=self.cfg.items_as('constants', float), srcex=self._src_exprs,
                       dt=dt)

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
            ufpts=self._scal_fpts_cpy, ILM=self.invLapMat, FLM=self.fluxLapMat
        )

    def makeSolLapMats(self, refmats):
        [Mxx, Mxy, Mxz, Myy, Myz, Mzz] = refmats
        rcpdjac = self.rcpdjac_at_np('upts') # (nupts, nelems)
        smats = self.smat_at_np('upts') # (ndims, nupts, ndims, nelems)

        [_, nelems] = np.shape(rcpdjac)
        invLapMat = np.zeros((self.nupts, self.nupts, nelems))

        if self.ndims == 2:
            for eidx in range(nelems):
                sxx = smats[0, :, 0, eidx]
                sxy = smats[1, :, 0, eidx]
                syx = smats[0, :, 1, eidx]
                syy = smats[1, :, 1, eidx]

                M =  (sxx*sxx + syx*syx)*Mxx # XX
                M += (sxx*sxy + syx*syy)*Mxy # XY
                M += (sxy*sxx + syy*syx)*Mxy # YX
                M += (sxy*sxy + syy*syy)*Myy # YY
                M = (M.T*rcpdjac[:, eidx]**2).T
                invLapMat[:,:,eidx] = np.linalg.inv(M)
        elif self.ndims == 3:
            for eidx in range(nelems):
                # Need to scale rows instead of columns? Use Transpose?
                # Smats correct order?
                sxx = smats[0, :, 0, eidx]
                sxy = smats[1, :, 0, eidx]
                sxz = smats[2, :, 0, eidx]
                syx = smats[0, :, 1, eidx]
                syy = smats[1, :, 1, eidx]
                syz = smats[2, :, 1, eidx]
                szx = smats[0, :, 2, eidx]
                szy = smats[1, :, 2, eidx]
                szz = smats[2, :, 2, eidx]

                M =  (sxx*sxx + syx*syx + szx*szx)*Mxx # XX
                M += (sxx*sxy + syx*syy + szx*szy)*Mxy # XY
                M += (sxx*sxz + syx*syz + szx*szz)*Mxz # XZ
                M += (sxy*sxx + syy*syx + szy*szx)*Mxy # YX
                M += (sxy*sxy + syy*syy + szy*szy)*Myy # YY
                M += (sxy*sxz + syy*syz + szy*szz)*Myz # YZ
                M += (sxz*sxx + syz*syx + szz*szx)*Mxz # ZX
                M += (sxz*sxy + syz*syy + szz*szy)*Myz # ZY
                M += (sxz*sxz + syz*syz + szz*szz)*Mzz # ZZ
                M = (M.T*rcpdjac[:, eidx]**2).T
                invLapMat[:,:,eidx] = np.linalg.inv(M)
        return invLapMat

    def makeFluxLapMats(self, refmats):
        [Mxx, Mxy, Mxz, Myy, Myz, Mzz] = refmats
        rcpdjac = self.rcpdjac_at_np('upts') # (nfpts, nelems)
        smats = self.smat_at_np('fpts') # (ndims, nfpts, ndims, nelems)

        [_, nelems] = np.shape(rcpdjac)
        LapMat = np.zeros((self.nupts, self.nfpts, nelems))

        if self.ndims == 2:
            for eidx in range(nelems):
                sxx = smats[0, :, 0, eidx]
                sxy = smats[1, :, 0, eidx]
                syx = smats[0, :, 1, eidx]
                syy = smats[1, :, 1, eidx]

                M =  (sxx*sxx + syx*syx)*Mxx # XX
                M += (sxx*sxy + syx*syy)*Mxy # XY
                M += (sxy*sxx + syy*syx)*Mxy # YX
                M += (sxy*sxy + syy*syy)*Myy # YY
                M = (M.T*rcpdjac[:, eidx]**2).T
                LapMat[:,:,eidx] = M
        elif self.ndims == 3:
            for eidx in range(nelems):
                # Need to scale rows instead of columns? Use Transpose?
                # Smats correct order?
                sxx = smats[0, :, 0, eidx]
                sxy = smats[1, :, 0, eidx]
                sxz = smats[2, :, 0, eidx]
                syx = smats[0, :, 1, eidx]
                syy = smats[1, :, 1, eidx]
                syz = smats[2, :, 1, eidx]
                szx = smats[0, :, 2, eidx]
                szy = smats[1, :, 2, eidx]
                szz = smats[2, :, 2, eidx]

                M =  (sxx*sxx + syx*syx + szx*szx)*Mxx # XX
                M += (sxx*sxy + syx*syy + szx*szy)*Mxy # XY
                M += (sxx*sxz + syx*syz + szx*szz)*Mxz # XZ
                M += (sxy*sxx + syy*syx + szy*szx)*Mxy # YX
                M += (sxy*sxy + syy*syy + szy*szy)*Myy # YY
                M += (sxy*sxz + syy*syz + szy*szz)*Myz # YZ
                M += (sxz*sxx + syz*syx + szz*szx)*Mxz # ZX
                M += (sxz*sxy + syz*syy + szz*szy)*Myz # ZY
                M += (sxz*sxz + syz*syz + szz*szz)*Mzz # ZZ
                M = (M.T*rcpdjac[:, eidx]**2).T
                LapMat[:,:,eidx] = M
        return LapMat
#  Interior Laplacian matrix
def refSolLapMat(self):
    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mxx = np.zeros((self.nupts, self.nupts)) 
    Mxy = np.zeros((self.nupts, self.nupts)) 
    Mxz = np.zeros((self.nupts, self.nupts)) if self.ndims == 3 else None
    Myy = np.zeros((self.nupts, self.nupts)) 
    Myz = np.zeros((self.nupts, self.nupts)) if self.ndims == 3 else None
    Mzz = np.zeros((self.nupts, self.nupts)) if self.ndims == 3 else None

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
                lagi = interpolate.lagrange(dfrpts, valsi)
                dlagi = lagi.deriv()
                ddlagi = lagi.deriv().deriv()

                valsj = np.zeros(p+3)
                valsj[j+1] = 1.0
                lagj = interpolate.lagrange(dfrpts, valsj)
                dlagj = lagj.deriv()
                ddlagj = lagj.deriv().deriv()


                for q in range(p+1):
                    for r in range(p+1):
                        Mxx[uidx(q,r), n] += ddlagi(upts_1d[q])*lagj(upts_1d[r])
                        Myy[uidx(q,r), n] += lagi(upts_1d[q])*ddlagj(upts_1d[r])
                        Mxy[uidx(q,r), n] += dlagi(upts_1d[q])*dlagj(upts_1d[r])

                n += 1


    elif self.ndims == 3:
        n = 0
        uidx = lambda i, j, k: i + j*(p+1) + k*(p+1)**2


        for k in range(p+1):
            for j in range(p+1):
                for i in range(p+1):
                    valsi = np.zeros(p+3)
                    valsi[i+1] = 1.0
                    lagi = interpolate.lagrange(dfrpts, valsi)
                    dlagi = lagi.deriv()
                    ddlagi = lagi.deriv().deriv()

                    valsj = np.zeros(p+3)
                    valsj[j+1] = 1.0
                    lagj = interpolate.lagrange(dfrpts, valsj)
                    dlagj = lagj.deriv()
                    ddlagj = lagj.deriv().deriv()

                    valsk = np.zeros(p+3)
                    valsk[k+1] = 1.0
                    lagk = interpolate.lagrange(dfrpts, valsk)
                    dlagk = lagk.deriv()
                    ddlagk = lagk.deriv().deriv()

                    for q in range(p+1):
                        for r in range(p+1):
                            for s in range(p+1):
                                Mxx[uidx(q,r,s), n] += ddlagi(upts_1d[q])*lagj(upts_1d[r])*lagk(upts_1d[s])
                                Mxy[uidx(q,r,s), n] += dlagi(upts_1d[q])*dlagj(upts_1d[r])*lagk(upts_1d[s])
                                Mxz[uidx(q,r,s), n] += dlagi(upts_1d[q])*lagj(upts_1d[r])*dlagk(upts_1d[s])
                                Myy[uidx(q,r,s), n] += lagi(upts_1d[q])*ddlagj(upts_1d[r])*lagk(upts_1d[s])
                                Myz[uidx(q,r,s), n] += lagi(upts_1d[q])*dlagj(upts_1d[r])*dlagk(upts_1d[s])
                                Mzz[uidx(q,r,s), n] += lagi(upts_1d[q])*lagj(upts_1d[r])*ddlagk(upts_1d[s])

                    n += 1

    return [Mxx, Mxy, Mxz, Myy, Myz, Mzz]

# Interface Laplacian matrix
def refFluxLapMat(self):
    p = self.basis.ubasis.order -1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack

    p = self.basis.ubasis.order-1
    upts_1d = self.basis.ubasis.pts[:p+1,0] # Hack
    Mxx = np.zeros((self.nupts, self.nfpts)) 
    Mxy = np.zeros((self.nupts, self.nfpts)) 
    Mxz = np.zeros((self.nupts, self.nfpts)) if self.ndims == 3 else None
    Myy = np.zeros((self.nupts, self.nfpts)) 
    Myz = np.zeros((self.nupts, self.nfpts)) if self.ndims == 3 else None
    Mzz = np.zeros((self.nupts, self.nfpts)) if self.ndims == 3 else None

    dfrpts = np.zeros(p+3)
    dfrpts[0] = -1.
    dfrpts[-1] = 1.
    dfrpts[1:-1] = upts_1d

    if self.ndims == 2:
        n = 0

        # QUAD FACE ORDERING:
        # 0 -> y = -1
        # 1 -> x =  1
        # 2 -> y =  1
        # 3 -> x = -1
        faceidx = lambda face, idx: face*(p+1) + idx
        uidx = lambda xidx, yidx: xidx + yidx*(p+1)

        for idx in range(p+1):
            valsa = np.zeros(p+3) # Negative faces
            valsa[0] = 1.0
            valsb = np.zeros(p+3) # Positive faces
            valsb[-1] = 1.0
            valsc = np.zeros(p+1) # Transverse
            valsc[idx] = 1.0

            laga = interpolate.lagrange(dfrpts, valsa)
            dlaga = laga.deriv()
            ddlaga = laga.deriv().deriv()
            lagb = interpolate.lagrange(dfrpts, valsb)
            dlagb = lagb.deriv()
            ddlagb = lagb.deriv().deriv()
            lagc = interpolate.lagrange(upts_1d, valsc)
            dlagc = lagc.deriv()
            ddlagc = lagc.deriv().deriv()

            for q in range(p+1):
                for r in range(p+1):
                    Mxx[uidx(q,r), faceidx(0, idx)] += ddlagc(upts_1d[q])*laga(upts_1d[r])
                    Mxx[uidx(q,r), faceidx(1, idx)] += ddlagb(upts_1d[q])*lagc(upts_1d[r])
                    Mxx[uidx(q,r), faceidx(2, idx)] += ddlagc(upts_1d[q])*lagb(upts_1d[r])
                    Mxx[uidx(q,r), faceidx(3, idx)] += ddlaga(upts_1d[q])*lagc(upts_1d[r])

                    Mxy[uidx(q,r), faceidx(0, idx)] += dlagc(upts_1d[q])*dlaga(upts_1d[r])
                    Mxy[uidx(q,r), faceidx(1, idx)] += dlagb(upts_1d[q])*dlagc(upts_1d[r])
                    Mxy[uidx(q,r), faceidx(2, idx)] += dlagc(upts_1d[q])*dlagb(upts_1d[r])
                    Mxy[uidx(q,r), faceidx(3, idx)] += dlaga(upts_1d[q])*dlagc(upts_1d[r])

                    Myy[uidx(q,r), faceidx(0, idx)] += lagc(upts_1d[q])*ddlaga(upts_1d[r])
                    Myy[uidx(q,r), faceidx(1, idx)] += lagb(upts_1d[q])*ddlagc(upts_1d[r])
                    Myy[uidx(q,r), faceidx(2, idx)] += lagc(upts_1d[q])*ddlagb(upts_1d[r])
                    Myy[uidx(q,r), faceidx(3, idx)] += laga(upts_1d[q])*ddlagc(upts_1d[r])


    elif self.ndims == 3:
        n = 0

        # HEX FACE ORDERING:
        # 0 -> z = -1
        # 1 -> y = -1
        # 2 -> x =  1
        # 3 -> y =  1
        # 4 -> x = -1
        # 5 -> z =  1
        faceidx = lambda face, a_idx, b_idx: face*(p+1)**2 + a_idx + b_idx*(p+1)
        uidx = lambda xidx, yidx, zidx: xidx + yidx*(p+1) + zidx*((p+1)**2)

        for aidx in range(p+1):
            for bidx in range(p+1):
                valsa = np.zeros(p+3) # Negative faces
                valsa[0] = 1.0
                valsb = np.zeros(p+3) # Positive faces
                valsb[-1] = 1.0
                valsc = np.zeros(p+1) # Transverse 1
                valsc[aidx] = 1.0
                valsd = np.zeros(p+1) # Transverse 2
                valsd[bidx] = 1.0

                laga = interpolate.lagrange(dfrpts, valsa)
                dlaga = laga.deriv()
                ddlaga = laga.deriv().deriv()

                lagb = interpolate.lagrange(dfrpts, valsb)
                dlagb = lagb.deriv()
                ddlagb = lagb.deriv().deriv()

                lagc = interpolate.lagrange(upts_1d, valsc)
                dlagc = lagc.deriv()
                ddlagc = lagc.deriv().deriv()

                lagd = interpolate.lagrange(upts_1d, valsd)
                dlagd = lagd.deriv()
                ddlagd = lagd.deriv().deriv()

                for q in range(p+1):
                    for r in range(p+1):
                        for s in range(p+1):
                            Mxx[uidx(q,r,s), faceidx(0, aidx, bidx)] += ddlagc(upts_1d[q])*lagd(upts_1d[r])*laga(upts_1d[s])
                            Mxx[uidx(q,r,s), faceidx(1, aidx, bidx)] += ddlagc(upts_1d[q])*laga(upts_1d[r])*lagd(upts_1d[s])
                            Mxx[uidx(q,r,s), faceidx(2, aidx, bidx)] += ddlagb(upts_1d[q])*lagc(upts_1d[r])*lagd(upts_1d[s])
                            Mxx[uidx(q,r,s), faceidx(3, aidx, bidx)] += ddlagc(upts_1d[q])*lagb(upts_1d[r])*lagd(upts_1d[s])
                            Mxx[uidx(q,r,s), faceidx(4, aidx, bidx)] += ddlaga(upts_1d[q])*lagc(upts_1d[r])*lagd(upts_1d[s])
                            Mxx[uidx(q,r,s), faceidx(5, aidx, bidx)] += ddlagc(upts_1d[q])*lagd(upts_1d[r])*lagb(upts_1d[s])
                    
                            Mxy[uidx(q,r,s), faceidx(0, aidx, bidx)] += dlagc(upts_1d[q])*dlagd(upts_1d[r])*laga(upts_1d[s])
                            Mxy[uidx(q,r,s), faceidx(1, aidx, bidx)] += dlagc(upts_1d[q])*dlaga(upts_1d[r])*lagd(upts_1d[s])
                            Mxy[uidx(q,r,s), faceidx(2, aidx, bidx)] += dlagb(upts_1d[q])*dlagc(upts_1d[r])*lagd(upts_1d[s])
                            Mxy[uidx(q,r,s), faceidx(3, aidx, bidx)] += dlagc(upts_1d[q])*dlagb(upts_1d[r])*lagd(upts_1d[s])
                            Mxy[uidx(q,r,s), faceidx(4, aidx, bidx)] += dlaga(upts_1d[q])*dlagc(upts_1d[r])*lagd(upts_1d[s])
                            Mxy[uidx(q,r,s), faceidx(5, aidx, bidx)] += dlagc(upts_1d[q])*dlagd(upts_1d[r])*lagb(upts_1d[s])
                            
                            Mxz[uidx(q,r,s), faceidx(0, aidx, bidx)] += dlagc(upts_1d[q])*lagd(upts_1d[r])*dlaga(upts_1d[s])
                            Mxz[uidx(q,r,s), faceidx(1, aidx, bidx)] += dlagc(upts_1d[q])*laga(upts_1d[r])*dlagd(upts_1d[s])
                            Mxz[uidx(q,r,s), faceidx(2, aidx, bidx)] += dlagb(upts_1d[q])*lagc(upts_1d[r])*dlagd(upts_1d[s])
                            Mxz[uidx(q,r,s), faceidx(3, aidx, bidx)] += dlagc(upts_1d[q])*lagb(upts_1d[r])*dlagd(upts_1d[s])
                            Mxz[uidx(q,r,s), faceidx(4, aidx, bidx)] += dlaga(upts_1d[q])*lagc(upts_1d[r])*dlagd(upts_1d[s])
                            Mxz[uidx(q,r,s), faceidx(5, aidx, bidx)] += dlagc(upts_1d[q])*lagd(upts_1d[r])*dlagb(upts_1d[s])
                            
                            Myy[uidx(q,r,s), faceidx(0, aidx, bidx)] += lagc(upts_1d[q])*ddlagd(upts_1d[r])*laga(upts_1d[s])
                            Myy[uidx(q,r,s), faceidx(1, aidx, bidx)] += lagc(upts_1d[q])*ddlaga(upts_1d[r])*lagd(upts_1d[s])
                            Myy[uidx(q,r,s), faceidx(2, aidx, bidx)] += lagb(upts_1d[q])*ddlagc(upts_1d[r])*lagd(upts_1d[s])
                            Myy[uidx(q,r,s), faceidx(3, aidx, bidx)] += lagc(upts_1d[q])*ddlagb(upts_1d[r])*lagd(upts_1d[s])
                            Myy[uidx(q,r,s), faceidx(4, aidx, bidx)] += laga(upts_1d[q])*ddlagc(upts_1d[r])*lagd(upts_1d[s])
                            Myy[uidx(q,r,s), faceidx(5, aidx, bidx)] += lagc(upts_1d[q])*ddlagd(upts_1d[r])*lagb(upts_1d[s])
                            
                            Myz[uidx(q,r,s), faceidx(0, aidx, bidx)] += lagc(upts_1d[q])*dlagd(upts_1d[r])*dlaga(upts_1d[s])
                            Myz[uidx(q,r,s), faceidx(1, aidx, bidx)] += lagc(upts_1d[q])*dlaga(upts_1d[r])*dlagc(upts_1d[s])
                            Myz[uidx(q,r,s), faceidx(2, aidx, bidx)] += lagb(upts_1d[q])*dlagc(upts_1d[r])*dlagc(upts_1d[s])
                            Myz[uidx(q,r,s), faceidx(3, aidx, bidx)] += lagc(upts_1d[q])*dlagb(upts_1d[r])*dlagc(upts_1d[s])
                            Myz[uidx(q,r,s), faceidx(4, aidx, bidx)] += laga(upts_1d[q])*dlagc(upts_1d[r])*dlagc(upts_1d[s])
                            Myz[uidx(q,r,s), faceidx(5, aidx, bidx)] += lagc(upts_1d[q])*dlagd(upts_1d[r])*dlagb(upts_1d[s])
                            
                            Mzz[uidx(q,r,s), faceidx(0, aidx, bidx)] += lagc(upts_1d[q])*lagd(upts_1d[r])*ddlaga(upts_1d[s])
                            Mzz[uidx(q,r,s), faceidx(1, aidx, bidx)] += lagc(upts_1d[q])*laga(upts_1d[r])*ddlagd(upts_1d[s])
                            Mzz[uidx(q,r,s), faceidx(2, aidx, bidx)] += lagb(upts_1d[q])*lagc(upts_1d[r])*ddlagd(upts_1d[s])
                            Mzz[uidx(q,r,s), faceidx(3, aidx, bidx)] += lagc(upts_1d[q])*lagb(upts_1d[r])*ddlagd(upts_1d[s])
                            Mzz[uidx(q,r,s), faceidx(4, aidx, bidx)] += laga(upts_1d[q])*lagc(upts_1d[r])*ddlagd(upts_1d[s])
                            Mzz[uidx(q,r,s), faceidx(5, aidx, bidx)] += lagc(upts_1d[q])*lagd(upts_1d[r])*ddlagb(upts_1d[s])
                        
    return [Mxx, Mxy, Mxz, Myy, Myz, Mzz]

