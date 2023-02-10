# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entropy'/>

<% inf = 1e20 %>
<%pyfr:kernel name='entropylocal' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='out fpdtype_t[${str(nfaces)}]'
              vb='in fpdtype_t[${str(nupts)}][2]'>
    // Compute minimum entropy across element
    fpdtype_t ui[${nvars}], vbi[2], d, p, e;

    fpdtype_t entmin = ${inf};
    for (int i = 0; i < ${nupts}; i++)
    {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        vbi[0] = vb[i][0];
        vbi[1] = vb[i][1];

        ${pyfr.expand('compute_entropy', 'ui', 'd', 'p', 'e', 'vbi')};

        entmin = fmin(entmin, e);
    }

    // Set interface entropy values to minimum
    % for i in range(nfaces):
    entmin_int[${i}] = entmin;
    % endfor
</%pyfr:kernel>
