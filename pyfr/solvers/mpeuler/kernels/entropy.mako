# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% inf = 1e20 %>
<%pyfr:macro name='compute_entropy' params='u, ad, amin, amax, p, e'>
    fpdtype_t d = ${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rcpd = 1.0/d, E = u[${ndims + nspec}];

    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${nspec + i}];
% endfor

    fpdtype_t a[${nspec}];
% for i in range(nspec - 1):
    a[${i}] = u[${nspec + ndims + 1 + i}];
% endfor
    a[${nspec - 1}] = 1 - (${' + '.join('a[{i}]'.format(i=i) for i in range(nspec - 1))});

    amin = ${inf};
    amax = ${-inf};
    ad = ${inf};
% for i in range(nspec):
    amin = min(amin, a[${i}]);
    amax = max(amax, a[${i}]);
    ad = min(ad, u[${i}]);
% endfor
    amax = 1 - amax;

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*rcpd*${pyfr.dot('rhov[{i}]', i=ndims)};
    fpdtype_t rcp_agm = 1/(${' + '.join('{rgm}*a[{i}]'.format(i=i, rgm=1/(c[f'gamma{i}']-1)) for i in range(nspec))});
    fpdtype_t gamma = rcp_agm + 1;
    p = rhoe*rcp_agm;

    // Compute combined species entropies
    e = ${' + '.join('u[{i}]*log(p*pow(a[{i}]/u[{i}], {gam}))'.format(i=i, gam=c[f'gamma{i}']) for i in range(nspec))};

    // Compute specific physical entropy
    e = ((p > 0) && (amin > 0) && (amax > 0) && (ad > 0)) ? exp(e) : ${inf};
</%pyfr:macro>
