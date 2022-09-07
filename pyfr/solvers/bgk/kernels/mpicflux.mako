# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='mpicflux' ndim='1'
              fl='inout view fpdtype_t[${str(nvars)}]'
              fr='in mpi fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'
              u='in broadcast fpdtype_t[${str(nvars)}][${str(ndims)}]'>
    // Perform the Riemann solve
    fpdtype_t Fn[${nvars}];
    ${pyfr.expand('rsolve', 'fl', 'fr', 'nl', 'Fn', 'u')};

    // Scale and write out the common normal fluxes
for (int i = 0; i < ${nvars}; i++) fl[i] = magnl*Fn[i];
</%pyfr:kernel>
