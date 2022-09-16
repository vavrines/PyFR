# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, q, f, v, p'>
    fpdtype_t phi = s[${ndims + 1}];
    p = s[0];
    fpdtype_t inv_d = 1/q[0];

    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = s[${i + 1}]*inv_d;
% endfor

% for i in range(ndims):
    // Mass flux
    f[${i}][0] = ${c['mpac-zeta']}*rhov[${i}];
    // Interface flux
    f[${i}][${ndims + 1}] = v[${i}]*phi;
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
