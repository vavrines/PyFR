# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='correct_pressure' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ucpy='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uinb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ILM='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              FLM='in fpdtype_t[${str(nupts)}][${str(nfpts)}]'>


    // Set divergence
    fpdtype_t divu[${nupts}];
    % for i in range(nupts):
        divu[${i}] = uoutb[${i}][${ndims}];
    % endfor

    // Compute pressure
    fpdtype_t q, P[${nupts}];

    // Calculate DFR pressure Poisson equation (negative interface contributions)
    % for i in range(nupts):
        q = (${' + '.join('ILM[{i}][{j}]*(divu[{j}])'.format(i=i, j=j) for j in range(nupts))})/${dt}; 
        P[${i}] = q + ucpy[${i}][${ndims}]; // Add to old pressure
    % endfor

    // Compute divF from forward Euler approximation of next state (store in uoutb)
    % for i,j in pyfr.ndrange(nupts, nvars-1):
        uoutb[${i}][${j}] = (uoutb[${i}][${j}] - uinb[${i}][${j}])/${dt};
    % endfor

    // Set pressure for uinb, set divergence to 0 for uoutb
    % for i in range(nupts):
        uinb[${i}][${nvars-1}] = P[${i}];
        uoutb[${i}][${nvars-1}] = 0.0;
    % endfor



</%pyfr:kernel>
