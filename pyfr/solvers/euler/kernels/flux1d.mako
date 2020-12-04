# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_1dflux' params='s, f, p, v'>
    fpdtype_t rho = s[0], invrho = 1.0/s[0], E = s[${nvars - 2}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

    // Compute the numerical entropy
    fpdtype_t ent = -rho*log(p*pow(rho, -${c['gamma']}))/(${c['gamma']- 1.0});
    ent = 0;

    // Density and energy fluxes
    f[0] = rhov[0];
    f[${nvars - 2}] = (E + p)*v[0];
    f[${nvars - 1}] = ent*v[0];

    // Momentum fluxes
% for i in range(ndims):
    f[${i + 1}] = rhov[0]*v[${i}]${' + p' if i == 0 else ''};
% endfor

</%pyfr:macro>
