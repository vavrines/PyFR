# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconffilter' ndim='1'
              t='scalar fpdtype_t'
              u_next='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'>

// Compute forward Euler approximation of -divF
% for i,v in pyfr.ndrange(nupts, nvars):
    u_next[${i}][${v}] = (u_next[${i}][${v}] - u[${i}][${v}])/${dt}; // Store in upts_outb
% endfor

</%pyfr:kernel>


