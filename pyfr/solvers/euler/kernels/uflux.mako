# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='uflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              f='out fpdtype_t[${str(ndims)}][${str(nvars)}]'>

// Compute the flux
fpdtype_t ftemp[${ndims}][${nvars}];

% for i, j in pyfr.ndrange(ndims, nvars):
    % if j == nvars - 2:
        ftemp[${i}][${j}] = ${c['gamma'] - 1}*(u[${j}] - 0.5*(pow(u[1], 2.0) + pow(u[2], 2.0))/u[0]); // Pressure
    % else:
        ftemp[${i}][${j}] = u[${j}];
    % endif
% endfor

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join('smats[{0}][{1}]*ftemp[{1}][{2}]'
                                 .format(i, k, j)
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
