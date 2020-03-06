# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<% gm = (c['gamma'] - 1.) %>
<% gp = (c['gamma'] + 1.) %>
<% rg = (1./c['gamma']) %>

// vn Leer Flux vector Splitting Approach
<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t fp[${nvars}],fm[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;
    
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Get the average, left, and right sound speeds
    fpdtype_t cl = sqrt(${c['gamma']}*pl/ul[0]);
    fpdtype_t cr = sqrt(${c['gamma']}*pr/ur[0]);

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Get normal mach numbers
    fpdtype_t ml = nvl/cl;
    fpdtype_t mr = nvr/cr;

    fpdtype_t ql = ${pyfr.dot('vl[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t qr = ${pyfr.dot('vr[{i}]', 'vr[{i}]', i=ndims)};


    // Get f+/- mass terms
    fpdtype_t fmp =  0.25*ul[0]*cl*(ml+1.)*(ml+1.);
    fpdtype_t fmm = -0.25*ur[0]*cr*(mr-1.)*(mr-1.);
    
    fp[0] = fmp;
    fm[0] = fmm;
% for i in range(ndims):
    fp[${i+1}] = fmp*(vl[${i}] + ${rg}*n[${i}]*(2*cl - nvl));
    fm[${i+1}] = fmm*(vr[${i}] - ${rg}*n[${i}]*(2*cr + nvr));
% endfor
    fp[${nvars-1}] = fmp*(0.5*(${gm}*nvl + 2*cl)*(${gm}*nvl + 2*cl)/${gm*gp} + 0.5*(ql - nvl*nvl));
    fm[${nvars-1}] = fmm*(0.5*(${gm}*nvr - 2*cr)*(${gm}*nvr - 2*cr)/${gm*gp} + 0.5*(qr - nvr*nvr));

    // Output
% for i in range(nvars):
    nf[${i}] = fm[${i}] + fp[${i}];
% endfor

</%pyfr:macro>