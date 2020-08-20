# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve3d' params='ul, ur, n, f'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    fpdtype_t a = sqrt(${0.25*c['gamma']}*(pl + pr)/(ul[0] + ur[0]))
                + 0.25*fabs(nv);

    // Output
% for i,j in pyfr.ndrange(ndims,nvars):
    f[${i}][${j}] = 0.5*(fl[${i}][${j}] + fr[${i}][${j}]) + n[${i}]*a*(ul[${j}] - ur[${j}]);
% endfor
</%pyfr:macro>
