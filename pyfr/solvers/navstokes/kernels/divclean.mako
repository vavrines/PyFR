# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='divclean' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              divu_unc='in fpdtype_t[${str(nupts)}]'
              divu_cor='in fpdtype_t[${str(nupts)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>


    // Calculate mesh spacing (assuming constant J: elements are cartesian)
    fpdtype_t h = 2.0*pow(1.0/rcpdjac[0], ${1.0/ndims});

    
    // Clean local divergence
    fpdtype_t U[${nupts*ndims}], RHS[${nupts}], Z[${nupts*ndims}];
    fpdtype_t u1mean[${ndims}] = {0}, u2mean[${ndims}] = {0};

    // Create velocity vector and compute cell-wise averages
    % for i,j in pyfr.ndrange(nupts, ndims):
        U[${i + j*nupts}] = uoutb[${i}][${j}];
        u1mean[${j}] += ${weights[i]}*(uoutb[${i}][${j}]);
    % endfor

    // Create RHS = -sum [ (du^I)*div(h(x_i))]
    % for i in range(nupts):
        RHS[${i}] = -(divu_cor[${i}] - divu_unc[${i}])*h;
    % endfor

    // Solve LSQ system
    % for i in range(nupts*ndims):
        Z[${i}] = (${' + '.join('{jx}*RHS[{j}]'.format(j=j, jx=jx)
                       for j, jx in enumerate(M[i]))});
    % endfor

    // Set velocity from 
    fpdtype_t du[${nupts*ndims}], tmp;
    % for i,j in pyfr.ndrange(nupts, ndims):
        tmp = Z[${i + j*nupts}] + 0*U[${i + j*nupts}];
        du[${i + j*nupts}] = (tmp - uoutb[${i}][${j}]);
        uoutb[${i}][${j}] = tmp;
    % endfor

    // Compute new cell-wise averages
    % for i,j in pyfr.ndrange(nupts, ndims):
        u2mean[${j}] += ${weights[i]}*(uoutb[${i}][${j}]);
    % endfor

    // Preserve means
    % for i,j in pyfr.ndrange(nupts, ndims):
        //uoutb[${i}][${j}] += u1mean[${j}] - u2mean[${j}];
    % endfor

    // Compute pressure correction
    fpdtype_t p1mean = 0, p2mean = 0;
    % for i in range(nupts):
        p1mean += ${weights[i]}*(uoutb[${i}][${ndims}]);
    % endfor

    fpdtype_t dp;
    % for i in range(nupts):
        dp = (${' + '.join('{jx}*du[{j}]'.format(j=j, jx=jx)
                       for j, jx in enumerate(P[i]))});
        //uoutb[${i}][${ndims}] -= dp;
    % endfor

    % for i in range(nupts):
        p2mean += ${weights[i]}*(uoutb[${i}][${ndims}]);
    % endfor

    % for i in range(nupts):
        //uoutb[${i}][${ndims}] += p1mean - p2mean;
    % endfor



</%pyfr:kernel>