# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='correct_pressure' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ufpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              ucpy='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uinb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ILM='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              FLM='in fpdtype_t[${str(nupts)}][${str(nfpts)}]'>


    // Set divergence
    // Calculate Laplacian of pressure from interface contributions
    fpdtype_t divu[${nupts}], intp[${nupts}];
    % for i in range(nupts):
        divu[${i}] = uoutb[${i}][${ndims}];
        intp[${i}] = (${' + '.join('FLM[{i}][{j}]*ufpts[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nfpts))})*${dt};
    % endfor

    // Compute pressure
    fpdtype_t q, P[${nupts}];

    // Calculate DFR pressure Poisson equation (negative interface contributions)
    % for i in range(nupts):
        q = (${' + '.join('ILM[{i}][{j}]*(divu[{j}] - intp[{j}])'.format(i=i, j=j) for j in range(nupts))})/${dt}; 
        P[${i}] = q; // Add to old pressure
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
