# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='limitinterp' ndim='2'
              u='inout fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              r='out fpdtype_t'>


    fpdtype_t E = u[${nvars-2}];

    u[0] = max(${tol}, u[0]);
    fpdtype_t invrho = 1.0/u[0];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
    % for i in range(ndims):
        rhov[${i}] = u[${i + 1}];
    % endfor

    fpdtype_t p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

    u[${nvars-2}] = p > ${tol} ? u[${nvars-2}] : ${tol/(c['gamma'] - 1)} + 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};

</%pyfr:kernel>

