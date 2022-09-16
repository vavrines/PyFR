# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mpaceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mpaceuler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ql='inout view fpdtype_t[${str(npass)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Compute the RHS
    fpdtype_t ur[${nvars}], qr[${npass}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'ql', 'norm_nl', 'ur', 'qr')};

    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ql', 'ur', 'qr' 'norm_nl', 'fn')};

    // Scale and write out the common normal fluxes
% for i in range(nvars):
    ul[${i}] = mag_nl*fn[${i}];
% endfor
</%pyfr:kernel>
