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
        outsoln = self.getAuxVariables(intg)

        solnfname    = self._writer.write(intg.soln, metadata, intg.tcurr)
        solnfnameaux = self._writeraux.write(outsoln, metadata, intg.tcurr)

        # If a post-action has been registered then invoke it
        self._invoke_postaction(mesh=intg.system.mesh.fname, soln=solnfname,
                                t=intg.tcurr)

        # Update the last output time
        self.tout_last = intg.tcurr

    def getAuxVariables(self, intg):
        elekeys = []
        outsoln = []
        for key in intg.system.ele_map:
            elekeys.append(key)
        for i in range(np.shape(intg.soln)[0]):
            soln = np.moveaxis(intg.soln[i],0,1)
            out = []
            f1 = intg.system.ele_map[elekeys[i]].F1._get()
            fk = intg.system.ele_map[elekeys[i]].fk._get()
            out.append(self.copyStructure(f1, soln[0]))
            out.append(self.copyStructure(fk, soln[0]))
            out = np.moveaxis(np.array(out),0,1)
            outsoln.append(out)
        return np.array(outsoln)

    def expandVarAcrossElem(self, var, targetshape):
        varshape = np.shape(var)
        out = np.zeros(targetshape)
        for i in range(targetshape[0]):
            out[i][:] = var
        return out

    def copyStructure(self, var, target):
        if np.shape(var) != np.shape(target):
            print('Varying shapes.', np.shape(var), np.shape(target))
        else:
            out = target*0
            for i in range(np.shape(var)[0]):
                for j in range(np.shape(var)[1]):
                        out[i][j] = var[i][j]
            return out