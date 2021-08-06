# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.kpp.kernels.flux'/>

<% tol = 1e-6 %>

<%pyfr:macro name='rsolve_f' params='ul, ur, n, nf'>
    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t nfl[${nvars}], nfr[${nvars}];
    fpdtype_t a, du;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    % for i in range(nvars):
        nfl[${i}] = (${' + '.join('n[{j}]*fl[{j}][{i}]'.format(i=i, j=j) for j in range(ndims))});
        nfr[${i}] = (${' + '.join('n[{j}]*fr[{j}][{i}]'.format(i=i, j=j) for j in range(ndims))});

        du  = ur[${i}] - ul[${i}];
        du  = abs(du) > ${tol} ? du : du > 0.0 ? ${tol} : ${-tol};
        a  = (nfr[${i}] - nfl[${i}])/du;
        nf[${i}] = a > 0 ? nfl[${i}] : nfr[${i}];
    % endfor

</%pyfr:macro>
<%pyfr:macro name='rsolve_b' params='ul, ur, n, nf'>
    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t nfl[${nvars}], nfr[${nvars}];
    fpdtype_t a, du;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    % for i in range(nvars):
        nfl[${i}] = (${' + '.join('n[{j}]*fl[{j}][{i}]'.format(i=i, j=j) for j in range(ndims))});
        nfr[${i}] = (${' + '.join('n[{j}]*fr[{j}][{i}]'.format(i=i, j=j) for j in range(ndims))});

        du  = ur[${i}] - ul[${i}];
        du  = abs(du) > ${tol} ? du : du > 0.0 ? ${tol} : ${-tol};
        a  = (nfr[${i}] - nfl[${i}])/du;
        nf[${i}] = a < 0 ? nfl[${i}] : nfr[${i}];
    % endfor

</%pyfr:macro>
