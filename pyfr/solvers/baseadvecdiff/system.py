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


        q1.enqueue(kernels['eles', 'enforce_positivity_interior'])
        # Copy solution to self._scal_upts_cpy
        q1.enqueue(kernels['eles', 'copy_soln'])

        # Take forward inviscid step
        q1.enqueue(kernels['eles', 'disu'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        q1.enqueue(kernels['eles', 'tdisf_inv'])
        q1.enqueue(kernels['eles', 'tdivtpcorf'])
        q1.enqueue(kernels['iint', 'comm_flux_inv_f'])
        q1.enqueue(kernels['bcint', 'comm_flux_inv_f'], t=t)
        q2.enqueue(kernels['mpiint', 'scal_fpts_send'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_unpack'])
        runall([q1, q2])
        q1.enqueue(kernels['mpiint', 'comm_flux_inv_f'])
        q1.enqueue(kernels['eles', 'tdivtconf'])
        q1.enqueue(kernels['eles', 'negdivconf_f'], t=t)

        # Take backward inviscid step
        q1.enqueue(kernels['eles', 'disu'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        q1.enqueue(kernels['eles', 'tdisf_inv'])
        q1.enqueue(kernels['eles', 'tdivtpcorf'])
        q1.enqueue(kernels['iint', 'comm_flux_inv_b'])
        q1.enqueue(kernels['bcint', 'comm_flux_inv_b'], t=t)
        q2.enqueue(kernels['mpiint', 'scal_fpts_send'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_unpack'])
        runall([q1, q2])
        q1.enqueue(kernels['mpiint', 'comm_flux_inv_b'])
        q1.enqueue(kernels['eles', 'tdivtconf'])
        q1.enqueue(kernels['eles', 'negdivconf_b'], t=t)

        # Take difference and switch with upts_cpy
        q1.enqueue(kernels['eles', 'get_du'])

        q1.enqueue(kernels['eles', 'disu'])
        q1.enqueue(kernels['eles', 'enforce_positivity_interface'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])

        if ('iint', 'copy_fpts') in kernels:
            q1.enqueue(kernels['iint', 'copy_fpts'])

        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
        if ('eles', 'shocksensor') in kernels:
            q1.enqueue(kernels['eles', 'shocksensor'])
            q1.enqueue(kernels['mpiint', 'artvisc_fpts_pack'])
        q1.enqueue(kernels['eles', 'tgradpcoru_upts'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_send'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_unpack'])

        runall([q1, q2])

        q1.enqueue(kernels['mpiint', 'con_u'])
        q1.enqueue(kernels['eles', 'tgradcoru_upts'])
        q1.enqueue(kernels['eles', 'gradcoru_upts'])
        q1.enqueue(kernels['eles', 'gradcoru_fpts'])

        q1.enqueue(kernels['mpiint', 'vect_fpts_pack'])
        if ('eles', 'shockvar') in kernels:
            q2.enqueue(kernels['mpiint', 'artvisc_fpts_send'])
            q2.enqueue(kernels['mpiint', 'artvisc_fpts_recv'])
            q2.enqueue(kernels['mpiint', 'artvisc_fpts_unpack'])

        runall([q1, q2])

        if ('eles', 'gradcoru_qpts') in kernels:
            q1.enqueue(kernels['eles', 'gradcoru_qpts'])
            
        q1.enqueue(kernels['eles', 'enforce_positivity_interior'])
        q1.enqueue(kernels['eles', 'tdisf_vis'])
        q1.enqueue(kernels['eles', 'tdivtpcorf'])
        q1.enqueue(kernels['iint', 'comm_flux_vis'])
        q1.enqueue(kernels['bcint', 'comm_flux_vis'], t=t)

        q2.enqueue(kernels['mpiint', 'vect_fpts_send'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_unpack'])

        runall([q1, q2])

        q1.enqueue(kernels['mpiint', 'comm_flux_vis'])
        q1.enqueue(kernels['eles', 'tdivtconf'])
        q1.enqueue(kernels['eles', 'negdivconf'], t=t)
        runall([q1])
