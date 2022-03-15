# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='clean_divergence' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    // Locally clean divergence
    fpdtype_t U[${nupts*ndims}], Ubar[${nupts*ndims}];
    % for i,j in pyfr.ndrange(nupts, ndims):
        U[${i + j*nupts}] = uoutb[${i}][${j}];
    % endfor
    % for i in range(nupts*ndims):
        Ubar[${i}] = (${' + '.join('{jx}*U[{j}]'.format(j=j, jx=jx)
                                 for j, jx in enumerate(V[i]) if jx != 0)});
    % endfor
    % for i,j in pyfr.ndrange(nupts, ndims):
        uoutb[${i}][${j}] = Ubar[${i + j*nupts}];
    % endfor




</%pyfr:kernel>
