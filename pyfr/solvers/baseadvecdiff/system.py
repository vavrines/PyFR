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
            
        self.compute_rhs(t, uinbank, foutbank)

    def compute_rhs(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active =   foutbank

        q1.enqueue(kernels['eles', 'disu'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])
        if ('eles', 'copy_soln') in kernels:
            q1.enqueue(kernels['eles', 'copy_soln'])
        if ('iint', 'copy_fpts') in kernels:
            q1.enqueue(kernels['iint', 'copy_fpts'])
        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
        runall([q1])
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
        runall([q1])

        q1.enqueue(kernels['eles', 'negdivconf_ns'], t=t)
        runall([q1])
        
        q1.enqueue(kernels['eles', 'tgradpcoru_unc'])
        q1.enqueue(kernels['eles', 'gradcoru_upts'])
        runall([q1])
        q1.enqueue(kernels['eles', 'compute_divergence_unc'])
        runall([q1])

        q1.enqueue(kernels['eles', 'disu_outb'])
        q1.enqueue(kernels['mpiint', 'scal_fpts_pack'])
        runall([q1])
        if ('eles', 'copy_soln') in kernels:
            q1.enqueue(kernels['eles', 'copy_soln'])
        if ('iint', 'copy_fpts') in kernels:
            q1.enqueue(kernels['iint', 'copy_fpts'])
        q1.enqueue(kernels['iint', 'con_u'])
        q1.enqueue(kernels['bcint', 'con_u'], t=t)
        q1.enqueue(kernels['eles', 'tgradpcoru_outb'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_send'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_recv'])
        q2.enqueue(kernels['mpiint', 'scal_fpts_unpack'])
        runall([q1, q2])
        q1.enqueue(kernels['mpiint', 'con_u'])
        q1.enqueue(kernels['eles', 'tgradcoru_upts'])
        q1.enqueue(kernels['eles', 'gradcoru_upts'])
        q1.enqueue(kernels['eles', 'gradcoru_fpts'])
        q1.enqueue(kernels['eles', 'compute_divergence_cor'])
        runall([q1, q2])

        
        
        q1.enqueue(kernels['eles', 'divclean'])
        q1.enqueue(kernels['eles', 'negdivconf_ns_inv'], t=t)
        runall([q1])
        
    def compute_grads(self, t, uinbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank

        q1.enqueue(kernels['eles', 'disu'])
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

        runall([q1])
