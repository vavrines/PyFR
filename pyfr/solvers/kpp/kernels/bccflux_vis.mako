# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.kpp.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.kpp.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux_vis' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    ${pyfr.expand('bc_common_flux_state_f', 'ul', 'gradul', 'artviscl', 'nl', 'magnl')};

</%pyfr:kernel>
