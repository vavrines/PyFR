# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
% for i, j in pyfr.ndrange(ndims, ndims):
    fout[${i}][${j}] += -${c['nu']}*grad_uin[${i}][${j}];
% endfor

% for i in range(ndims):
    fout[${i}][${ndims}] = 0.0;
% endfor
</%pyfr:macro>
