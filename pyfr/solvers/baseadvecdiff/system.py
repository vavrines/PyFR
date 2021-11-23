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

        piter = self.cfg.getint('solver', 'piter', 1)
        for _ in range(piter):
            self.iterate_pressure(t, uinbank, foutbank)
            
        self.compute_rhs(t, uinbank, foutbank)

    def iterate_pressure(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        # Get interface common pressures
        q1.enqueue(kernels['eles', 'copy_soln'])
        q1.enqueue(kernels['eles', 'disu'])
        runall([q1])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])
        if ('iint', 'copy_fpts') in kernels:
            q1.enqueue(kernels['iint', 'copy_fpts'])
        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
        q1.enqueue(kernels['eles', 'tgradpcoru_upts'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_send'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_unpack'])
        runall([q1, q2])
        q1.enqueue(kernels['mpiint', 'con_u'])
        q1.enqueue(kernels['eles', 'copy_fpts2'])


        # Compute intermediate state (without pressure component)
        q1.enqueue(kernels['eles', 'zero_interior_pressure'])
        q1.enqueue(kernels['eles', 'disu'])
        runall([q1])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])
        if ('iint', 'copy_fpts') in kernels:
            q1.enqueue(kernels['iint', 'copy_fpts'])
        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
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

        runall([q1, q2])

        if ('eles', 'gradcoru_qpts') in kernels:
            q1.enqueue(kernels['eles', 'gradcoru_qpts'])
        q1.enqueue(kernels['eles', 'tdisf'])
        q1.enqueue(kernels['eles', 'tdivtpcorf'])
        q1.enqueue(kernels['iint', 'comm_flux'])
        q1.enqueue(kernels['bcint', 'comm_flux'], t=t)

        q2.enqueue(kernels['mpiint', 'vect_fpts_send'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_unpack'])

        runall([q1, q2])

        q1.enqueue(kernels['mpiint', 'comm_flux'])
        q1.enqueue(kernels['eles', 'tdivtconf'])
        q1.enqueue(kernels['eles', 'negdivconf_ns'], t=t)
        runall([q1])


        # Compute new divergence -------------
        q1.enqueue(kernels['eles', 'disu_outb'])
        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
        if ('eles', 'shocksensor') in kernels:
            q1.enqueue(kernels['eles', 'shocksensor'])
            q1.enqueue(kernels['mpiint', 'artvisc_fpts_pack'])
        q1.enqueue(kernels['eles', 'tgradpcoru_upts_outb'])
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

        q1.enqueue(kernels['eles', 'compute_divergence'])
        runall([q1])
        # -------------------

        q1.enqueue(kernels['eles', 'correct_pressure'])
        runall([q1])

    def compute_rhs(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        q1.enqueue(kernels['eles', 'disu'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])

        if ('eles', 'copy_soln') in kernels:
            q1.enqueue(kernels['eles', 'copy_soln'])
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
        q1.enqueue(kernels['eles', 'tdisf'])
        q1.enqueue(kernels['eles', 'tdivtpcorf'])
        q1.enqueue(kernels['iint', 'comm_flux'])
        q1.enqueue(kernels['bcint', 'comm_flux'], t=t)

        q2.enqueue(kernels['mpiint', 'vect_fpts_send'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'vect_fpts_unpack'])

        runall([q1, q2])

        q1.enqueue(kernels['mpiint', 'comm_flux'])
        q1.enqueue(kernels['eles', 'tdivtconf'])
        q1.enqueue(kernels['eles', 'negdivconf'], t=t)
        runall([q1])