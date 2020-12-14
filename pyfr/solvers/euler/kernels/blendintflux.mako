# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='blendintflux' ndim='1'
              f_HO='inout fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              f_LO='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              alpha_fpts='in fpdtype_t[${str(nfpts)}]'
              >


% for i, var in pyfr.ndrange(nfpts, nvars):
    f_HO[${i}][${var}] = alpha_fpts[${i}]*f_LO[${i}][${var}] + (1.0 - alpha_fpts[${i}])*f_HO[${i}][${var}];
% endfor


</%pyfr:kernel>
