# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='divclean' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    
    // Clean local divergence
    fpdtype_t U[${nupts*ndims}], RHS[${nupts}], Ubar[${nupts*ndims}];
    fpdtype_t u1mean[${ndims}] = {0}, u2mean[${ndims}] = {0};

    % for i,j in pyfr.ndrange(nupts, ndims):
        U[${i + j*nupts}] = uoutb[${i}][${j}];
        u1mean[${j}] += ${weights[i]}*(uoutb[${i}][${j}]);
    % endfor

    % for i in range(nupts):
        RHS[${i}] = -(${' + '.join('{jx}*(U[{j}])'.format(j=j, jx=jx)
                                 for j, jx in enumerate(D[i]) if jx != 0)});
    % endfor

    % for i in range(nupts*ndims):
        Ubar[${i}] = (${' + '.join('{jx}*(RHS[{j}])'.format(j=j, jx=jx)
                       for j, jx in enumerate(M[i]))});
    % endfor

    % for i,j in pyfr.ndrange(nupts, ndims):
        uoutb[${i}][${j}] = Ubar[${i + j*nupts}] + U[${i + j*nupts}];
    % endfor

    % for i,j in pyfr.ndrange(nupts, ndims):
        u2mean[${j}] += ${weights[i]}*(uoutb[${i}][${j}]);
    % endfor

    % for i,j in pyfr.ndrange(nupts, ndims):
        uoutb[${i}][${j}] += u1mean[${j}] - u2mean[${j}];
    % endfor
</%pyfr:kernel>