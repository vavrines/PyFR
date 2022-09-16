# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpaceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ql, ur, qr, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}], pl, pr;

    ${pyfr.expand('inviscid_flux', 'ul', 'ql', 'fl', 'vl', 'pl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'qr', 'fr', 'vr', 'pr')};

    // Normal of the average interface velocity
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Estimate the wave speed
    fpdtype_t al = fabs(nvl) + sqrt(nvl*nvl + ql[0]*${c['mpac-zeta']});
    fpdtype_t ar = fabs(nvr) + sqrt(nvr*nvr + qr[0]*${c['mpac-zeta']});
    fpdtype_t a  = max(al, ar);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
