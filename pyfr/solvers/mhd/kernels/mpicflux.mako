# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>

<% beta, tau = c['ldg-beta'], c['ldg-tau'] %>

<%pyfr:kernel name='mpicflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout mpi fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              gradur='in mpi fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              artviscr='in mpi fpdtype_t'
              entminl='in view fpdtype_t'
              entminr='in mpi fpdtype_t'
              entmin_intl='inout view fpdtype_t'
              entmin_intr='inout mpi fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor

entmin_intl = entmin_intr = fmin(entminl, entminr);
</%pyfr:kernel>
