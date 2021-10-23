# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'ul[{i}] + ur[{i}]', i=ndims)};

    if (nv >= 0) {
        % for i in range(nvars):
            nf[${i}] = (${' + '.join('n[{j}]*fl[{j}][{i}]'
                                 .format(i=i, j=j) for j in range(ndims))});
        % endfor
    }
    else {
        % for i in range(nvars):
            nf[${i}] = (${' + '.join('n[{j}]*fr[{j}][{i}]'
                                 .format(i=i, j=j) for j in range(ndims))});
        % endfor
    }
	

</%pyfr:macro>
