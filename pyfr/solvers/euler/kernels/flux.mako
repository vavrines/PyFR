<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, p, v, ploc'>
    fpdtype_t invrho = 1.0/s[0], E = s[${nvars - 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];
    f[${i}][${nvars - 1}] = (E + p)*v[${i}];
% endfor

    // Momentum fluxes
    fpdtype_t vb[${ndims}] = {0};
    fpdtype_t r = sqrt(ploc[0]*ploc[0] + ploc[1]*ploc[1]);
    vb[0] = -${c['omg']}*r*ploc[1];
    vb[1] =  ${c['omg']}*r*ploc[0];
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*(v[${j}] - vb[${j}])${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
