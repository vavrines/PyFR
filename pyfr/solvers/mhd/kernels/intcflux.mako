# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>

<% beta, tau = c['ldg-beta'], c['ldg-tau'] %>

<%pyfr:kernel name='intcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              gradur='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              artviscr='in view fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] =  magnl*(ficomm[${i}]);
    ur[${i}] = -magnl*(ficomm[${i}]);
% endfor

artviscl = fmin(artviscl, artviscr);
artviscr = fmin(artviscl, artviscr);
</%pyfr:kernel>
