# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='enforce_positivity' ndim='2'
              u='inout fpdtype_t[${str(nvars)}]'>

<% tol = 1e-6 %>
fpdtype_t E = u[${nvars-1}];

u[0] = fmax(${tol}, u[0]);
fpdtype_t invrho = 1.0/u[0];

// Compute the velocities
fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${i + 1}];
% endfor

fpdtype_t p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

u[${nvars-1}] = p > ${tol} ? u[${nvars-1}] : ${tol/(c['gamma'] - 1)} + 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};

</%pyfr:kernel>

