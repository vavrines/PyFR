# -*- coding: utf-8 -*-
<%inherit file='base'/>

<%pyfr:macro name='entropy' params='u, s'>
    fpdtype_t rho = u[0], invrho = 1.0/u[0], E = u[${nvars - 2}];

    fpdtype_t rhov[${ndims}];
    % for i in range(ndims):
        rhov[${i}] = u[${i + 1}];
    % endfor

    // Compute the pressure
    fpdtype_t p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

    // Compute the numerical entropy
    s = -rho*log(p*pow(rho, -${c['gamma']}))/(${c['gamma']- 1.0});
</%pyfr:macro>

<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='residual' ndim='2'
              tdivtconf='in fpdtype_t[${str(nvars)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              r='out fpdtype_t'>

    fpdtype_t u_next[${nvars}], s0, s1;

    % for i in range(nvars):
        u_next[${i}] = u[${i}] - ${dt}*rcpdjac*tdivtconf[${i}];
    % endfor

    ${pyfr.expand('entropy', 'u', 's0')};
    ${pyfr.expand('entropy', 'u_next', 's1')};

    fpdtype_t dsdt = (s1 - s0)/${dt};
    fpdtype_t dfdx = rcpdjac*tdivtconf[${nvars-1}];

    r = dsdt + dfdx;

</%pyfr:kernel>

