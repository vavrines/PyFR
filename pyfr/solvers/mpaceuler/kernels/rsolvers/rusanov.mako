# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpaceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ql, ur, qr, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'ql', 'fl', 'vl', 'pl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'qr', 'fr', 'vr', 'pr')};

    fpdtype_t vl[${ndims}] = ${pyfr.array('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t vr[${ndims}] = ${pyfr.array('ur[{i}]', i=(1, ndims + 1))};

    // Normal of the average interface velocity
    fpdtype_t nv = 0.5*${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the wave speed
    fpdtype_t a = fabs(nv) + sqrt(nv*nv + ${c['ac-zeta']});

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
