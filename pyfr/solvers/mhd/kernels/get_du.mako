# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='get_du' ndim='2'
              u_orig='inout fpdtype_t[${str(nvars)}]'
              u_rev='inout fpdtype_t[${str(nvars)}]'>

fpdtype_t du;

% for i in range(nvars):
	du = u_rev[${i}] - u_orig[${i}];
    u_rev[${i}] = u_orig[${i}];
    u_orig[${i}] = du;
% endfor

</%pyfr:kernel>
