# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux1d'/>

// E-CUSP Approxiamte Riemann Solver
<% al0 = 0.05 %>
<%pyfr:macro name='rsolve_t1d' params='ul, ur, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_1dflux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_1dflux', 'ur', 'fr', 'pr', 'vr')};

    // Compute Roe averaged states
    //fpdtype_t roa = sqrt(ul[0])*sqrt(ur[0]);
    //fpdtype_t ha = (sqrt(ul[0])*(pr + ur[${ndims + 1}]) + sqrt(ur[0])*(pl + ul[${ndims + 1}]))
    //                    /(sqrt(ul[0])*ur[0] + sqrt(ur[0])*ul[0]);
    //fpdtype_t invsqrulpur = 1/(sqrt(ul[0]) + sqrt(ur[0]));

% for i in range(ndims):
    ##va[${i}] = (vl[${i}]*sqrt(ul[0]) + vr[${i}]*sqrt(ur[0]))*invsqrulpur;
% endfor

    //fpdtype_t qq = ${pyfr.dot('va[{i}]', 'va[{i}]', i=ndims)};
    //fpdtype_t a = sqrt(${c['gamma'] - 1}*(ha - 0.5*qq));
    // Estimate the maximum wave speed
    fpdtype_t a = sqrt(${c['gamma']}*(pl + pr)/(ul[0] + ur[0]));
    fpdtype_t ma = 0.5*(vl[0] + vr[0])/a;
    fpdtype_t ama = fabs(ma);

    fpdtype_t be = (ama >= 1.) ?  ma/ama : (ma >= 0.) ? max(0.,2.*ma - 1.) :
                                                        min(0.,2.*ma + 1.);

    fpdtype_t al = (ama < ${al0}) ? 0.5*(${al0} + ama*ama*${1./al0}) : ama;
    fpdtype_t as = al*a - 0.5*be*(vl[0] + vr[0]);

    //printf("be: %E \n",be);

% for i in range(nvars):
    nf[${i}] = 0.5*((1.+be)*fl[${i}] + (1.-be)*fr[${i}]) - 0.5*as*(ur[${i}] - ul[${i}]);
% endfor

</%pyfr:macro>

<%include file='pyfr.solvers.euler.kernels.rsolvers.transformed'/>
