from collections import defaultdict

import numpy as np

from pyfr.mpiutil import get_comm_rank_root, mpi
from pyfr.plugins.base import BasePlugin, init_csv


class YPlusPlugin(BasePlugin):
    name = 'yplus'
    systems = ['navier-stokes']
    formulations = ['dual', 'std']

    def __init__(self, intg, cfgsect, suffix):
        super().__init__(intg, cfgsect, suffix)

        comm, rank, root = get_comm_rank_root()

        # Output frequency
        self.nsteps = self.cfg.getint(cfgsect, 'nsteps')

        # Get bins for y+ histogram
        self.bins = self.cfg.getliteral(cfgsect, 'bins') + [np.inf]
        if 0 not in self.bins:
            self.bins = [0] + self.bins
        self.bins = np.array(self.bins, dtype=float)
        if not np.all(np.diff(self.bins) > 0):
            raise ValueError ('Bins for y+ histogram must be monotonically'
                              'increasing.')

        self.nbins = len(self.bins) - 1
        self.normalize_bins = self.cfg.getbool(cfgsect, 'normalize-bins',
                                               False)

        # Viscous correction
        self._viscorr = self.cfg.get('solver', 'viscosity-correction', 'none')

        # Constant variables
        self._constants = self.cfg.items_as('constants', float)

        # Underlying elements class
        self.elementscls = intg.system.elementscls

        # Boundary to integrate over
        bc = f'bcon_{suffix}_p{intg.rallocs.prank}'

        # Get the mesh and elements
        mesh, elemap = intg.system.mesh, intg.system.ele_map

        # See which ranks have the boundary
        bcranks = comm.gather(bc in mesh, root=root)

        # The root rank needs to open the output file
        if rank == root:
            if not any(bcranks):
                raise RuntimeError(f'Boundary {suffix} does not exist')

            # CSV header
            header = ['t', 'yp_min', 'yp_max']
            header += [f'<{b}' for b in self.bins[1:-1]]
            header += [f'>{self.bins[-2]}']
            header += ['x_min', 'y_min', 'z_min'][:self.ndims]
            header += ['x_max', 'y_max', 'z_max'][:self.ndims]

            # Open
            self.outf = init_csv(self.cfg, cfgsect, ','.join(header))

        # Interpolation matrices and quadrature weights
        self._m0 = m0 = {}
        self._m4 = m4 = {}
        rcpjact = {}
        h = {}
        ploc_fpts = {}
        self._qwts = qwts = defaultdict(list)

        # If we have the boundary then process the interface
        if bc in mesh:
            # Element indices, associated face normals and relative flux
            # points position with respect to the moments origin
            eidxs = defaultdict(list)
            norms = defaultdict(list)

            for etype, eidx, fidx, flags in mesh[bc].astype('U4,i4,i1,i2'):
                eles = elemap[etype]
                
                # Phyiscal normals
                pnorms = eles.get_pnorms(eidx, fidx)
                pnnorms = pnorms/np.linalg.norm(pnorms, axis=1)[:, np.newaxis]

                if (etype, fidx) not in m0:
                    facefpts = eles.basis.facefpts[fidx]

                    m0[etype, fidx] = eles.basis.m0[facefpts]
                    qwts[etype, fidx] = eles.basis.fpts_wts[facefpts]

                if etype not in m4:
                    m4[etype] = eles.basis.m4

                    # Get the smats at the solution points
                    smat = eles.smat_at_np('upts').transpose(2, 0, 1, 3)

                    # Get |J|^-1 at the solution points
                    rcpdjac = eles.rcpdjac_at_np('upts')

                    # Product to give J^-T at the solution points
                    rcpjact[etype] = smat*rcpdjac

                    # Get element-wise average flux point locations at boundary
                    pfpts = eles.ploc_at_np('fpts')[facefpts].T
                    ploc_fpts[etype] = np.dot(pfpts, qwts[etype, fidx]).T
                    ploc_fpts[etype] /= np.sum(qwts[etype, fidx])

                    # Get wall normal mesh spacing
                    fsmats = eles.smat_at_np('fpts')[:,facefpts]
                    h[etype] = 2*np.linalg.norm(np.einsum('jkil,kj->ikl', fsmats, pnnorms), axis=0)

                eidxs[etype, fidx].append(eidx)
                norms[etype, fidx].append(pnorms)

            self._eidxs = {k: np.array(v) for k, v in eidxs.items()}
            self._norms = {k: np.array(v) for k, v in norms.items()}

            self._rcpjact = {k: rcpjact[k[0]][..., v]
                                for k, v in self._eidxs.items()}
            self._h = {k: h[k[0]][..., v]
                          for k, v in self._eidxs.items()}
            self._ploc_fpts = {k: ploc_fpts[k[0]][..., v]
                                  for k, v in self._eidxs.items()}

    def __call__(self, intg):
        # Return if no output is due
        if intg.nacptsteps % self.nsteps:
            return

        # MPI info
        comm, rank, root = get_comm_rank_root()

        # Solution matrices indexed by element type
        solns = dict(zip(intg.system.ele_types, intg.soln))
        ndims, nvars = self.ndims, self.nvars

        mu = self._constants['mu']
        yp_min = np.inf
        yp_max = 0.0
        ploc_yp_min, ploc_yp_max = np.zeros((2, self.ndims))

        hist = np.zeros(self.nbins, dtype=int)

        for etype, fidx in self._m0:
            # Get the interpolation operator
            m0 = self._m0[etype, fidx]
            nfpts, nupts = m0.shape

            # Extract the relevant elements from the solution
            uupts = solns[etype][..., self._eidxs[etype, fidx]]

            # Interpolate to the face
            ufpts = m0 @ uupts.reshape(nupts, -1)
            ufpts = ufpts.reshape(nfpts, nvars, -1)
            ufpts = ufpts.swapaxes(0, 1)

            # Get the quadrature weights and normal vectors
            qwts = self._qwts[etype, fidx]
            norms = self._norms[etype, fidx]
            
            # Get operator and J^-T matrix
            m4 = self._m4[etype]
            rcpjact = self._rcpjact[etype, fidx]

            # Transformed gradient at solution points
            tduupts = m4 @ uupts.reshape(nupts, -1)
            tduupts = tduupts.reshape(ndims, nupts, nvars, -1)

            # Physical gradient at solution points
            duupts = np.einsum('ijkl,jkml->ikml', rcpjact, tduupts)
            duupts = duupts.reshape(ndims, nupts, -1)

            # Interpolate gradient to flux points
            dufpts = np.array([m0 @ du for du in duupts])
            dufpts = dufpts.reshape(ndims, nfpts, nvars, -1)
            dufpts = dufpts.swapaxes(1, 2)

            # Get velocity gradients
            gradrho, gradrhou = dufpts[:, 0], dufpts[:, 1:-1]
            rho = ufpts[0]
            gradu = (gradrhou - gradrho[:, None]*ufpts[None, 1:-1]/rho) / rho

            # Compute wall shear stress
            tau_wall =  mu*np.linalg.norm(np.einsum('ijkl,lki->jkl', gradu, norms), axis=0)

            # Compute friction velocity
            u_tau = np.sqrt(tau_wall/rho)

            # Get wall distance
            h = self._h[etype, fidx]

            # Compute y+
            y_plus = u_tau*h*rho/mu

            # Average y+ across face
            y_plus = np.dot(y_plus.T, qwts)

            # Bin y+ values
            hist += np.histogram(y_plus, self.bins)[0]

            # Find min/max y+
            minidx = np.argmin(y_plus)
            if y_plus[minidx] < yp_min:
                yp_min = y_plus[minidx]
                ploc_yp_min = self._ploc_fpts[etype, fidx][:, minidx]
                
            maxidx = np.argmax(y_plus)
            if y_plus[maxidx] > yp_max:
                yp_max = y_plus[maxidx]
                ploc_yp_max = self._ploc_fpts[etype, fidx][:, maxidx]
        
        # Reduce
        if rank != root:
            comm.reduce((yp_min, ploc_yp_min), op=mpi.MIN, root=root)
            comm.reduce((yp_max, ploc_yp_max), op=mpi.MAX, root=root)
            comm.Reduce(hist, None, op=mpi.SUM, root=root)
        else:
            yp_min, ploc_yp_min = comm.reduce((yp_min, ploc_yp_min), op=mpi.MIN, root=root)
            yp_max, ploc_yp_max = comm.reduce((yp_max, ploc_yp_max), op=mpi.MAX, root=root)
            comm.Reduce(mpi.IN_PLACE, hist, op=mpi.SUM, root=root)

            # Normalize bin counts
            if self.normalize_bins:
                hist = hist/np.sum(hist)

            fm = [yp_min, yp_max] + list(hist) + list(ploc_yp_min) + list(ploc_yp_max)

            # Write
            print(intg.tcurr, *fm, sep=',', file=self.outf)

            # Flush to disk
            self.outf.flush()
