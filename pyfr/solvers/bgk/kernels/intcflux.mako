# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='intcflux' ndim='1'
              fl='inout view fpdtype_t[${str(nvars)}]'
              fr='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    // Perform the Riemann solve
    fpdtype_t Fn[${nvars}];
    ${pyfr.expand('rsolve', 'fl', 'fr', 'nl', 'Fn')};

    // Scale and write out the common normal fluxes
% for i in range(nvars):
    fl[${i}] =  magnl*Fn[${i}];
    fr[${i}] = -magnl*Fn[${i}];
% endfor
</%pyfr:kernel>
