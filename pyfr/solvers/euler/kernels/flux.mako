# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f'>
    fpdtype_t p = s[${ndims}];

    // Momentum fluxes
    % for i, j in pyfr.ndrange(ndims, ndims):
        f[${i}][${j}] = s[${i}]*s[${j}]${' + p' if i == j else ''};
    % endfor

    % for i in range(ndims):
        f[${i}][${ndims}] = 0.0;
    % endfor
</%pyfr:macro>
