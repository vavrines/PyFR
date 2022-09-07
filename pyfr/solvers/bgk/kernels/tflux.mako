# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='tflux' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nvars)}][${str(ndims)}]'>
    // Compute and transform the fluxes
    fpdtype_t ftemp[${ndims}];
for (int j = 0; j < ${nvars}; j++) {
    % for i in range(ndims):
    ftemp[${i}] = -u[j][${i}]*f[j]; // Flux = -u(x,y,t).f(x,y,t)
    % endfor

    % for i in range(ndims):
    F[${i}][j] = ${' + '.join(f'smats[{i}][{k}]*ftemp[{k}]' for k in range(ndims))};
    % endfor
}

</%pyfr:kernel>
