# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;
    fpdtype_t rl = 1.0/ul[0], rr = 1.0/ur[0];
    fpdtype_t nf_fl, nf_fr, nf_flstar, nf_frstar, d_star;
    fpdtype_t rcp_lstar, rcp_rstar;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Compute the Roe-averaged velocity
    fpdtype_t nv = (sqrt(rl)*nvl + sqrt(rr)*nvr)
                 / (sqrt(rl) + sqrt(rr));

    // Compute the Roe-averaged enthalpy
    fpdtype_t H = (sqrt(rl)*(pr + ur[${ndims + 1}])
                 + sqrt(rr)*(pl + ul[${ndims + 1}]))
                / (sqrt(rl)*rr + sqrt(rr)*rl);

    // Roe average sound speed
    fpdtype_t a = sqrt(${c['gamma'] - 1}*(H - 0.5*nv*nv));

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = nv - a;
    fpdtype_t sr = nv + a;
    fpdtype_t s_star = (pr - pl + rl*nvl*(sl - nvl) -
                        rr*nvr*(sr - nvr)) /
                       (rl*(sl - nvl) - rr*(sr - nvr));

    // Output
    rcp_lstar = 1 / (sl - s_star);
    rcp_rstar = 1 / (sr - s_star);
% for i in range(nvars):
    nf_fl = ${' + '.join('n[{j}]*fl[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
    nf_fr = ${' + '.join('n[{j}]*fr[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
% if i == 0:
    nf_flstar = s_star*(sl*ul[${i}] - nf_fl) * rcp_lstar;
    nf_frstar = s_star*(sr*ur[${i}] - nf_fr) * rcp_rstar;
% else:
    d_star = ${'s_star' if i == nvars - 1 else 'n[{0}]'.format(i - 1)};
    nf_flstar = (s_star*(sl*ul[${i}] - nf_fl) +
                 sl*(pl + rl*(sl - nvl)*(s_star - nvl))*d_star) *
                rcp_lstar;
    nf_frstar = (s_star*(sr*ur[${i}] - nf_fr) +
                 sr*(pr + rr*(sr - nvr)*(s_star - nvr))*d_star) *
                rcp_rstar;
% endif

    nf[${i}] = (sl >= 0) ? nf_fl : (sl <= 0 && s_star >= 0) ? nf_flstar :
               (s_star <= 0 && sr >= 0) ? nf_frstar : nf_fr;
% endfor
</%pyfr:macro>
