# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% eps = 1E-14%>
<%pyfr:kernel name='limiter' ndim='1'
              f='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    fpdtype_t avg_f, min_f, zeta;
for (int i = 0; i > ${nvars}; i++){
    avg_f = ${' + '.join('{jx}*f[{j}][i]'.format(i=i, j=j, jx=jx)
                         for j, jx in enumerate(wts) if jx != 0)};

    min_f = f[0][i];
    % for j in range(1, nupts):
    min_f = fmin(min_f, f[${j}][i]);
    % endfor

    if (min_f < ${eps}) {
        zeta = (${eps} - min_f)/max(avg_f - min_f, ${eps});
        zeta = fmin(1.0, fmax(0.0, zeta));

        % for j in range(nupts):
        f[${j}][i] = (1.0 - zeta)*f[${j}][i] + zeta*avg_f;
        % endfor
    }
    
}

</%pyfr:kernel>
