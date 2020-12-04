# -*- coding: utf-8 -*-

from pyfr.inifile import Inifile
from pyfr.plugins.base import BasePlugin
from pyfr.writers.native import NativeWriter
import numpy as np


class WriterPlugin(BasePlugin):
    name = 'writer'
    systems = ['*']
    formulations = ['dual', 'std']

    def __init__(self, intg, cfgsect, suffix=None):
        super().__init__(intg, cfgsect, suffix)

        # Construct the solution writer
        basedir = self.cfg.getpath(cfgsect, 'basedir', '.', abs=True)
        basename = self.cfg.get(cfgsect, 'basename')
        self._writer = NativeWriter(intg, self.nvars, basedir, basename,
                                    prefix='soln')
        self._writeraux = NativeWriter(intg, 2, basedir, basename,
                                    prefix='aux')

        # Output time step and last output time
        self.dt_out = self.cfg.getfloat(cfgsect, 'dt-out')
        self.tout_last = intg.tcurr

        # Output field names
        self.fields = intg.system.elementscls.convarmap[self.ndims]

        # Register our output times with the integrator
        intg.call_plugin_dt(self.dt_out)

        # If we're not restarting then write out the initial solution
        if not intg.isrestart:
            self.tout_last -= self.dt_out
            self(intg)

    def __call__(self, intg):
        if intg.tcurr - self.tout_last < self.dt_out - self.tol:
            return

        stats = Inifile()
        stats.set('data', 'fields', ','.join(self.fields))
        stats.set('data', 'prefix', 'soln')
        intg.collect_stats(stats)

        # Prepare the metadata
        metadata = dict(intg.cfgmeta,
                        stats=stats.tostr(),
                        mesh_uuid=intg.mesh_uuid)

        # Write out the file
        solnfname = self._writer.write(intg.soln, metadata, intg.tcurr)
        nsects, nupts, nvars, nelems = shape = np.shape(intg.soln)
        outsoln = self.getAuxVariables(intg, shape)
        solnfnameaux = self._writeraux.write(outsoln, metadata, intg.tcurr)

        # If a post-action has been registered then invoke it
        self._invoke_postaction(mesh=intg.system.mesh.fname, soln=solnfname,
                                t=intg.tcurr)

        # Update the last output time
        self.tout_last = intg.tcurr

    def getAuxVariables(self, intg, shape):
        nsects, nupts, nvars, nelems = shape
        elekeys = []
        outsoln = []
        for key in intg.system.ele_map:
            elekeys.append(key)
        for i in range(nsects):
            r = intg.system.ele_map[elekeys[i]].residual._get()
            r = r[:,np.newaxis, :]
            outsoln.append(r)
        return np.array(outsoln)


