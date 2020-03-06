# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux1d'/>

<%pyfr:macro name='rsolve_t1d' params='ul, ur, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_1dflux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_1dflux', 'ur', 'fr', 'pr', 'vr')};

    // Estimate the maximum wave speed
    fpdtype_t s = sqrt(${c['gamma']}*(pl + pr)/(ul[0] + ur[0]))
                + 0.5*fabs(vl[0] + vr[0]);

    //fpdtype_t al = sqrt(${c['gamma']}*pl/ul[0]);
    //fpdtype_t ar = sqrt(${c['gamma']}*pr/ur[0]);
    //fpdtype_t s = max(fabs(vl[0])+al,fabs(vr[0])+ar);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(fl[${i}] + fr[${i}]) - 0.5*s*(ur[${i}] - ul[${i}]);
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.euler.kernels.rsolvers.transformed'/>