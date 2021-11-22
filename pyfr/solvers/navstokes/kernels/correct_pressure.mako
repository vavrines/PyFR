# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='correct_pressure' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uinb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ucpy='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ufpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              ILM='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              FLM='in fpdtype_t[${str(nupts)}][${str(nfpts)}]'>

    // Set RHS and calculate Laplacian of pressure from interface contributions
    fpdtype_t RHS[${nupts}], intp[${nupts}];
    % for i in range(nupts):
        RHS[${i}] = uoutb[${i}][${ndims}]/${dt} + ucpy[${i}][1] - ucpy[${i}][0]; // divu + Lap(p^n) - Lap(nu*divu)
        intp[${i}] = (${' + '.join('FLM[{i}][{j}]*ufpts[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nfpts))});
    % endfor


    // Calculate DFR pressure Poisson equation (negative interface contributions)
    % for i in range(nupts):
        uinb[${i}][${nvars-1}] = (${' + '.join('ILM[{i}][{j}]*(RHS[{j}] - intp[{j}])'.format(i=i, j=j) for j in range(nupts))});
        //uinb[${i}][${nvars-1}] = (${' + '.join('ILM[{i}][{j}]*(RHS[{j}])'.format(i=i, j=j) for j in range(nupts))});
    % endfor


</%pyfr:kernel>
