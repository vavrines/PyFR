# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='add_visc' ndim='2'
              grads='inout fpdtype_t[${ndims}][${str(nvars)}]'
              revvisc='in fpdtype_t[${str(nvars)}]'>

% for i, j in pyfr.ndrange(ndims, nvars):
	grads[${i}][${j}] *= revvisc[${j}];
% endfor

</%pyfr:kernel>
