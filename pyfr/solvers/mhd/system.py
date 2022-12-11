from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.mhd.elements import MHDElements
from pyfr.solvers.mhd.inters import (MHDIntInters, MHDMPIInters,
                                     MHDBaseBCInters)
from pyfr.util import memoize


class MHDSystem(BaseAdvectionSystem):
    name = 'mhd'

    elementscls = MHDElements
    intinterscls = MHDIntInters
    mpiinterscls = MHDMPIInters
    bbcinterscls = MHDBaseBCInters

    @memoize
    def _rhs_graphs(self, uinbank, foutbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, foutbank)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # # Interpolate the solution to the flux points
        g1.add_all(k['eles/disu'])

        # Pack and send these interpolated solutions to our neighbours
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disu'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # If entropy filtering, pack and send the entropy values to neighbors
        g1.add_mpi_reqs(m['ent_fpts_recv'])
        g1.add_all(k['mpiint/ent_fpts_pack'])
        for send, pack in zip(m['ent_fpts_send'], k['mpiint/ent_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute common entropy minima at internal/boundary interfaces
        g1.add_all(k['iint/comm_entropy'],
                   deps=k['mpiint/ent_fpts_pack'])
        g1.add_all(k['bcint/comm_entropy'],
                   deps=k['eles/disu'])

        # Compute the common normal flux at our internal/boundary interfaces
        g1.add_all(k['iint/comm_flux'],
                   deps=k['eles/disu'] + k['mpiint/scal_fpts_pack'])
        g1.add_all(k['bcint/comm_flux'],
                   deps=k['eles/disu'] + k['bcint/comm_entropy'])

        # Make a copy of the solution (if used by source terms)
        g1.add_all(k['eles/copy_soln'])

        # Interpolate the solution to the quadrature points
        g1.add_all(k['eles/qptsu'])

        # Compute the transformed flux
        for l in k['eles/tdisf_curved'] + k['eles/tdisf_linear']:
            ldeps = deps(l, 'eles/qptsu')
            g1.add(l, deps=ldeps)

        # Compute the transformed divergence of the partially corrected flux
        for l in k['eles/tdivtpcorf']:
            ldeps = deps(l, 'eles/tdisf_curved', 'eles/tdisf_linear',
                         'eles/copy_soln', 'eles/disu')
            g1.add(l, deps=ldeps + k['mpiint/scal_fpts_pack'])
        g1.commit()

        g2 = self.backend.graph()

        # Compute the common normal flux at our MPI interfaces
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/comm_flux']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Compute common entropy minima at MPI interfaces
        g2.add_all(k['mpiint/ent_fpts_unpack'])
        for l in k['mpiint/comm_entropy']:
            g2.add(l, deps=deps(l, 'mpiint/ent_fpts_unpack'))

        # Compute the transformed divergence of the corrected flux
        g2.add_all(k['eles/tdivtconf'], deps=k['mpiint/comm_flux'])

        # Obtain the physical divergence of the corrected flux
        for l in k['eles/negdivconf']:
            g2.add(l, deps=deps(l, 'eles/tdivtconf'))
        g2.commit()

        return g1, g2

    @memoize
    def _preproc_graphs(self, uinbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, None)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()

        # Interpolate the solution to the flux points
        if 'eles/local_entropy' in k:
            g1.add_all(k['eles/disu'])

        # Compute local minimum entropy within element
        g1.add_all(k['eles/local_entropy'])

        # Pack and send the entropy values to neighbors
        g1.add_all(k['mpiint/ent_fpts_pack'], deps=k['eles/local_entropy'])
        for send, pack in zip(m['ent_fpts_send'], k['mpiint/ent_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute common entropy minima at internal/boundary interfaces
        g1.add_all(k['iint/comm_entropy'], deps=k['eles/local_entropy'])
        g1.add_all(k['bcint/comm_entropy'],
                   deps=k['eles/local_entropy'] + k['eles/disu'])
        g1.commit()

        if 'mpiint/comm_entropy' in k:
            # Compute common entropy minima at MPI interfaces
            g2 = self.backend.graph()

            g2.add_all(k['mpiint/ent_fpts_unpack'])
            for l in k['mpiint/comm_entropy']:
                g2.add(l, deps=deps(l, 'mpiint/ent_fpts_unpack'))
            g2.commit()

            return g1, g2
        else:
            return g1,

    def postproc(self, uinbank):
        k, _ = self._get_kernels(uinbank, None)

        self.backend.run_kernels(k['eles/entropy_filter'])

    @memoize
    def _powell_graphs(self, uinbank, foutbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, foutbank)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # # Interpolate the solution to the flux points
        g1.add_all(k['eles/disu'])

        # Pack and send these interpolated solutions to our neighbours
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disu'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute the common normal flux at our internal/boundary interfaces
        g1.add_all(k['iint/comm_flux'],
                   deps=k['eles/disu'] + k['mpiint/scal_fpts_pack'])
        g1.add_all(k['bcint/comm_flux'],
                   deps=k['eles/disu'])

        # Make a copy of the solution (if used by source terms)
        g1.add_all(k['eles/copy_soln'])

        # Interpolate the solution to the quadrature points
        g1.add_all(k['eles/qptsu'])

        # Compute the transformed flux
        for l in k['eles/tdisf_curved'] + k['eles/tdisf_linear']:
            ldeps = deps(l, 'eles/qptsu')
            g1.add(l, deps=ldeps)

        # Compute the transformed divergence of the partially corrected flux
        for l in k['eles/tdivtpcorf']:
            ldeps = deps(l, 'eles/tdisf_curved', 'eles/tdisf_linear',
                         'eles/copy_soln', 'eles/disu')
            g1.add(l, deps=ldeps + k['mpiint/scal_fpts_pack'])
        g1.commit()

        g2 = self.backend.graph()

        # Compute the common normal flux at our MPI interfaces
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/comm_flux']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Compute the transformed divergence of the corrected flux
        g2.add_all(k['eles/tdivtconf'], deps=k['mpiint/comm_flux'])

        # Obtain the physical divergence of the corrected flux
        for l in k['eles/powellsource']:
            g2.add(l, deps=deps(l, 'eles/tdivtconf'))
        g2.commit()

        return g1, g2
