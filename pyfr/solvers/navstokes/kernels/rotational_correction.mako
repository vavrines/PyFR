# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='rotational_correction' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ucpy='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ufpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              SLM='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              FLM='in fpdtype_t[${str(nupts)}][${str(nfpts)}]'>


    // Calculate Laplacian of divu from interface/solution contributions
    fpdtype_t ldu[${nupts}];
    % for i in range(nupts):
        ldu[${i}] =  (${' + '.join('FLM[{i}][{j}]*ufpts[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nfpts))});
        ldu[${i}] += (${' + '.join('SLM[{i}][{j}]*uoutb[{j}][{var}]'.format(i=i, j=j, var=ndims) for j in range(nupts))});
    % endfor

    // Store in ucpy[0]
    % for i in range(nupts):
        ucpy[${i}][0] = ${c['nu']}*ldu[${i}]; 
    % endfor


</%pyfr:kernel>
