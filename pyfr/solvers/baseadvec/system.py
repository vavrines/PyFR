# -*- coding: utf-8 -*-

from pyfr.solvers.base import BaseSystem


class BaseAdvectionSystem(BaseSystem):
    _nqueues = 2

    def rhs(self, t, uinbank, foutbank):
        self.rhs_HO(t, uinbank, foutbank)

    def rhs_LO(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        # Calculate low-order gradients
        q1 << kernels['eles', 'disu_LO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_LO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('eles', 'copy_soln_at_fpts') in kernels:
            q1 << kernels['eles', 'copy_soln_at_fpts']()

        q1 << kernels['eles', 'tdisf']()
        q1 << kernels['eles', 'tdivtpcorf_LO']()

        q1 << kernels['iint', 'comm_flux_LO']()
        q1 << kernels['bcint', 'comm_flux_LO'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_LO']()

        q1 << kernels['eles', 'tdivtconf_LO']()
        q1 << kernels['eles', 'residual']()


        # Map u to f shape and calculate gradient
        q1 << kernels['eles', 'udisf']()
        q1 << kernels['eles', 'tdivtpcorf_LO']()

        q1 << kernels['iint', 'comm_flux_LO']()
        q1 << kernels['bcint', 'comm_flux_LO'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_LO']()
        q1 << kernels['eles', 'tdivtconf_LO']()



        # q1 << kernels['eles', 'normalizeresidual']()


        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_HO_int']()
        q1 << kernels['eles', 'limitinterp']()

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('eles', 'copy_soln_at_fpts') in kernels:
            q1 << kernels['eles', 'copy_soln_at_fpts']()
        q1 << kernels['iint', 'comm_flux']()
        q1 << kernels['bcint', 'comm_flux'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux']()

        q1 << kernels['eles', 'riemanndifference']()
        q1 << kernels['eles', 'tdivtpcorf_RD']()
        q1 << kernels['eles', 'tdivtconf_RD']()

        if ('eles', 'tdivf_qpts') in kernels:
            raise NotImplementedError()
        else:
            q1 << kernels['eles', 'negdivconf'](t=t)
        runall([q1])

    def rhs_HO(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])

        q1 << kernels['eles', 'disu_HO_int']()
        q1 << kernels['eles', 'limitinterp']()

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('eles', 'copy_soln_at_fpts') in kernels:
            q1 << kernels['eles', 'copy_soln_at_fpts']()
        q1 << kernels['eles', 'tdisf']()

        q1 << kernels['eles', 'tdivtpcorf_HO']()

        q1 << kernels['iint', 'comm_flux_LO']()
        q1 << kernels['bcint', 'comm_flux_LO'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()

        runall([q1, q2])

        q1 << kernels['mpiint', 'comm_flux_LO']()

        q1 << kernels['eles', 'tdivtconf_HO']()
        q1 << kernels['eles', 'residual']()
        q1 << kernels['eles', 'normalizeresidual']()



        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_HO_int']()
        q1 << kernels['eles', 'limitinterp']()

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('eles', 'copy_soln_at_fpts') in kernels:
            q1 << kernels['eles', 'copy_soln_at_fpts']()
        q1 << kernels['iint', 'comm_flux']()
        q1 << kernels['bcint', 'comm_flux'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux']()

        q1 << kernels['eles', 'riemanndifference']()
        q1 << kernels['eles', 'tdivtpcorf_RD']()
        q1 << kernels['eles', 'tdivtconf_RD']()

        if ('eles', 'tdivf_qpts') in kernels:
            raise NotImplementedError()
        else:
            q1 << kernels['eles', 'negdivconf'](t=t)
        runall([q1])
