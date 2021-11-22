# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='old_pressure_laplacian' ndim='1'
              ucpy='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ufpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              SLM='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              FLM='in fpdtype_t[${str(nupts)}][${str(nfpts)}]'>


    // Calculate Laplacian of old pressure from interface/solution contributions
    fpdtype_t lp[${nupts}];
    % for i in range(nupts):
        lp[${i}] =  (${' + '.join('FLM[{i}][{j}]*ufpts[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nfpts))});
        lp[${i}] += (${' + '.join('SLM[{i}][{j}]*ucpy[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nupts))});
    % endfor

    // Store in ucpy[1]
    % for i in range(nupts):
        ucpy[${i}][1] = lp[${i}]; 
    % endfor


</%pyfr:kernel>
