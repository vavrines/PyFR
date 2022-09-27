# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, a, d, p, v'>
    d = ${' + '.join('s[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t invrho = 1/d, E = s[${ndims + nspec}];

% for i in range(nspec - 1):
    a[${i}] = s[${nspec + ndims + 1 + i}];
% endfor
    a[${nspec - 1}] = 1 - (${' + '.join('a[{i}]'.format(i=i) for i in range(nspec - 1))});

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + nspec}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};
    fpdtype_t agm = ${' + '.join('{rgm}*a[{i}]'.format(i=i, rgm=1/(c[f'gamma{i}']-1)) for i in range(nspec))};
    p = rhoe/agm;

    // Mass flux
% for i, j in pyfr.ndrange(ndims, nspec):
    f[${i}][${j}] = s[${j}]*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + nspec}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor

    // Energy fluxes
% for i in range(ndims):
    f[${i}][${ndims + nspec}] = (E + p)*v[${i}];
% endfor

    // Species flux (zeroed as set in negdivconf)
% for i, j in pyfr.ndrange(ndims, nspec-1):
    f[${i}][${j + nspec + ndims + 1}] = 0.;
% endfor

</%pyfr:macro>
