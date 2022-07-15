# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entropy'/>

<% inf = 1e20 %>
<% niters = 20 %>

<%pyfr:macro name='get_minima' params='u, dmin, pmin, emin'>
    fpdtype_t d, p, e;
    fpdtype_t ui[${nvars}];

    dmin = ${inf}, pmin = ${inf}, emin = ${inf};

    for (int i = 0; i < ${nupts}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        ${pyfr.expand('compute_entropy', 'ui', 'd', 'p', 'e')};
        dmin = fmin(dmin, d); pmin = fmin(pmin, p); emin = fmin(emin, e);
    }
    % if viscous:
    emin = ${inf}; // Do not apply entropy constraints on viscous component
    % endif
</%pyfr:macro>

<%pyfr:macro name='apply_filter' params='umodes, ffac, uf, zeta'>
    // Precompute filter factors
    % for i in range(nupts):
    ffac[${i}] = exp(-zeta*${ubdegs2[i]});
    % endfor

    // Compute filtered solution
    for (int uidx = 0; uidx < ${nupts}; uidx++) {
        % for j in range(nvars):
        uf[uidx][${j}] = 0.0;
        % endfor

        for (int midx = 0; midx < ${nupts}; midx++) {
            for (int vidx = 0; vidx < ${nvars}; vidx++) {
                tmp = ffac[midx]*umodes[midx][vidx]; // Filtered mode
                uf[uidx][vidx] += vdm[uidx][midx]*tmp;
            }
        }
    }
</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin='in fpdtype_t'
              vdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'>
    fpdtype_t dmin, pmin, emin;

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'dmin', 'pmin', 'emin')};

    // Filter if out of bounds
    if ((dmin < ${d_min}) || (pmin < ${p_min}) || (emin < entmin - ${e_tol})) {
        // Compute modal basis
        fpdtype_t umodes[${nupts}][${nvars}] = {{0}};

        for (int uidx = 0; uidx < ${nupts}; uidx++) {
            for (int vidx = 0; vidx < ${nvars}; vidx++) {
                for (int midx = 0; midx < ${nupts}; midx++) {
                    umodes[uidx][vidx] += invvdm[uidx][midx]*u[midx][vidx];
                }
            }
        }

        // Setup filter
        fpdtype_t zeta_low = 0.0;
        fpdtype_t zeta_high = ${zeta_max};
        fpdtype_t tmp, zeta;

        fpdtype_t uf[${nupts}][${nvars}] = {{0}};
        fpdtype_t ffac[${nupts}];

        // Iterate filter strength with bisection algorithm
        for (int iter = 0; iter < ${niters}; iter++) {
            // Bisection estimate
            zeta = 0.5*(zeta_low + zeta_high);

            ${pyfr.expand('apply_filter', 'umodes', 'ffac', 'uf', 'zeta')};
            ${pyfr.expand('get_minima', 'uf', 'dmin', 'pmin', 'emin')};

            // Compute new bracket
            if ((dmin < ${d_min}) || (pmin < ${p_min}) || (emin < entmin - ${e_tol})) {
                zeta_low = zeta;
            }
            else {
                zeta_high = zeta;
            }
        }

        // Apply filtered solution with bounds-preserving filter strength
        zeta = zeta_high;
        ${pyfr.expand('apply_filter', 'umodes', 'ffac', 'u', 'zeta')};
    }
    
</%pyfr:kernel>
