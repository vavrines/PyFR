# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='density' ndim='2'
              u='inout fpdtype_t[${str(nvars)}]',
              q='inout fpdtype_t[${str(npass)}]'>
    fpdtype_t phi = u[${ndims + 1}];
    fpdtype_t d_new = ${0.5*(c['rho0'] + c['rho1'])} - ${0.5*(c['rho0'] - c['rho1'])}*tanh(${c['mpac-sigma']}*phi);
    fpdtype_t dnrdo = d_new/q[0];

% for i in range(ndims):
    u[${i + 1}] = u[${i + 1}]*dnrdo;
% endfor
    q[0] = d_new;
</%pyfr:kernel>