<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, p, v, vb'>
    fpdtype_t invrho = 1.0/s[0], E = s[${nvars - 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    fpdtype_t rote = 0.5*s[0]*(vb[0]*vb[0] + vb[1]*vb[1]);
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)} + rote);

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];
    f[${i}][${nvars - 1}] = (E + p)*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${j}]*v[${i}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
