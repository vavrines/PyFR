# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr;
    fpdtype_t al[${nspec}], ar[${nspec}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'al', 'rhol', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'ar', 'rhor', 'pr', 'vr')};

    // phase averaged gamma
    fpdtype_t gl = (1/${' + '.join('al[{i}]*{rgm}').format(i=i, rgm=1/(c[f'gamma{i}']-1)) for i in range(nspec)}) + 1;
    fpdtype_t gr = (1/${' + '.join('ar[{i}]*{rgm}').format(i=i, rgm=1/(c[f'gamma{i}']-1)) for i in range(nspec)}) + 1;

    // Sum the left and right velocities and take the normal
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // wave speeds
    fpdtype_t cl = sqrt(gl*pl/rhol);
    fpdtype_t cr = sqrt(gr*pr/rhor);
    fpdtype_t c = min(fabs(nvl) + cl, fasb(nvr) + cr);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + c*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
