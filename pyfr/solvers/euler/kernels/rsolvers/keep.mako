# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux1d'/>

<%pyfr:macro name='v_and_p' params='u,v,p'>
    fpdtype_t invrho = 1./u[0];

% for i in range(ndims):
    v[${i}] = invrho*u[${i+1}];
% endfor
    p = ${c['gamma'] - 1}*(u[${nvars-1}] - 0.5*u[0]*${pyfr.dot('v[{i}]', i=ndims)});

</%pyfr:macro>

<% rgm = 1./(c['gamma'] - 1) %>
<%pyfr:macro name='rsolve_t1d' params='ul, ur, nf'>
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('v_and_p','ul','vl','pl')};
    ${pyfr.expand('v_and_p','ur','vr','pr')};

    fpdtype_t rb = 0.5*(ul[0] + ur[0]);
    fpdtype_t ub = 0.5*(vl[0] + vr[0]);

    nf[0] = rb*ub;
% for i in range(ndims):
    nf[${i+1}] = 0.5*rb*ub*(vl[${i}] + vr[${i}]);
% endfor
   nf[${nvars-1}] = 0.5*ub*(rb*(${pyfr.dot('vl[{i}]', 'vr[{i}]', i=ndims)} +
                                    ${rgm}*(pl/ul[0] + pr/ur[0])) + (pl + pr));

</%pyfr:macro>
