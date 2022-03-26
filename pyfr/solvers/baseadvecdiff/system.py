# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem


class BaseAdvectionDiffusionSystem(BaseAdvectionSystem):
    def rhs(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        q1 << kernels['eles', 'disu_ext']()
        q1 << kernels['eles', 'calcentropy']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])

        q1 << kernels['eles', 'disu_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('iint', 'copy_fpts') in kernels:
            q1 << kernels['iint', 'copy_fpts']()
        q1 << kernels['iint', 'con_u']()
        q1 << kernels['bcint', 'con_u'](t=t)
        if ('eles', 'shocksensor') in kernels:
            q1 << kernels['eles', 'shocksensor']()
            q1 << kernels['mpiint', 'artvisc_fpts_pack']()
        q1 << kernels['mpiint', 'entmin_fpts_pack']()
        q1 << kernels['mpiint', 'entminint_fpts_pack']()

        q1 << kernels['eles', 'tgradpcoru_upts']()
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()

        runall([q1, q2])

        q1 << kernels['mpiint', 'con_u']()
        q1 << kernels['eles', 'tgradcoru_upts_ext']()
        q1 << kernels['eles', 'gradcoru_upts_ext']()
        q1 << kernels['eles', 'gradcoru_fpts_ext']()
        q1 << kernels['mpiint', 'vect_fpts_pack']()
        if ('eles', 'shockvar') in kernels:
            q2 << kernels['mpiint', 'artvisc_fpts_send']()
            q2 << kernels['mpiint', 'artvisc_fpts_recv']()
            q2 << kernels['mpiint', 'artvisc_fpts_unpack']()

        q2 << kernels['mpiint', 'entmin_fpts_send']()
        q2 << kernels['mpiint', 'entmin_fpts_recv']()
        q2 << kernels['mpiint', 'entmin_fpts_unpack']()
        q2 << kernels['mpiint', 'entminint_fpts_send']()
        q2 << kernels['mpiint', 'entminint_fpts_recv']()
        q2 << kernels['mpiint', 'entminint_fpts_unpack']()
        runall([q1, q2])

        q1 << kernels['eles', 'tgradcoru_upts_int']()
        q1 << kernels['eles', 'gradcoru_upts_int']()
        q1 << kernels['eles', 'gradcoru_fpts_int']()
        if ('eles', 'gradcoru_qpts') in kernels:
            q1 << kernels['eles', 'gradcoru_qpts']()

        # Split viscous and inviscid parts
        # Viscous
        q1 << kernels['eles', 'tdisf_vis']()
        q1 << kernels['eles', 'tdivtpcorf']()
        q1 << kernels['iint', 'comm_flux_vis']()
        q1 << kernels['bcint', 'comm_flux_vis'](t=t)
        q2 << kernels['mpiint', 'vect_fpts_send']()
        q2 << kernels['mpiint', 'vect_fpts_recv']()
        q2 << kernels['mpiint', 'vect_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_vis']()
        q1 << kernels['eles', 'tdivtconf']()
        runall([q1])
        q1 << kernels['eles', 'copy_divf']()

        # Inviscid
        q1 << kernels['eles', 'disu_ext']()
        q1 << kernels['eles', 'disu_int']()
        q1 << kernels['mpiint', 'vect_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'tdisf_inv']()
        q1 << kernels['eles', 'tdivtpcorf']()
        q1 << kernels['iint', 'comm_flux_inv']()
        q1 << kernels['bcint', 'comm_flux_inv'](t=t)
        q2 << kernels['mpiint', 'vect_fpts_send']()
        q2 << kernels['mpiint', 'vect_fpts_recv']()
        q2 << kernels['mpiint', 'vect_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_inv']()
        q1 << kernels['eles', 'tdivtconf']()
        runall([q1])

        # Filter output
        q1 << kernels['eles', 'calcminentropy']()
        runall([q1])

        q1 << kernels['eles', 'negdivconf'](t=t)
        runall([q1])

    def compute_grads(self, t, uinbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank

        q1 << kernels['eles', 'disu_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_int']()
        q1 << kernels['iint', 'con_u']()
        q1 << kernels['bcint', 'con_u'](t=t)

        q1 << kernels['eles', 'tgradpcoru_upts']()
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()

        runall([q1, q2])

        q1 << kernels['mpiint', 'con_u']()
        q1 << kernels['eles', 'tgradcoru_upts_ext']()
        q1 << kernels['eles', 'gradcoru_upts_ext']()
        q1 << kernels['eles', 'gradcoru_fpts_ext']()

        q1 << kernels['eles', 'tgradcoru_upts_int']()
        q1 << kernels['eles', 'gradcoru_upts_int']()
        q1 << kernels['eles', 'gradcoru_fpts_int']()



        # q1 << kernels['eles', 'disu']()
        # q1 << kernels['mpiint', 'scal_fpts_pack']()
        # runall([q1])

        # if ('iint', 'copy_fpts') in kernels:
        #     q1 << kernels['iint', 'copy_fpts']()
        # q1 << kernels['iint', 'con_u']()
        # q1 << kernels['bcint', 'con_u'](t=t)
        # q1 << kernels['eles', 'tgradpcoru_upts']()
        # q2 << kernels['mpiint', 'scal_fpts_send']()
        # q2 << kernels['mpiint', 'scal_fpts_recv']()
        # q2 << kernels['mpiint', 'scal_fpts_unpack']()

        # runall([q1, q2])

        # q1 << kernels['mpiint', 'con_u']()
        # q1 << kernels['eles', 'tgradcoru_upts_ext']()
        # q1 << kernels['eles', 'gradcoru_upts_ext']()
        # runall([q1])
        # q1 << kernels['eles', 'tgradcoru_upts_int']()
        # q1 << kernels['eles', 'gradcoru_upts_int']()

        runall([q1])