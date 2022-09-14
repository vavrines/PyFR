# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%pyfr:kernel name='limiter' ndim='1'
              f='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>
<% eps = 1e-12%>

fpdtype_t avg_f, min_f, beta;
for (int i = 0; i > ${nvars}; i++){
    avg_f = ${' + '.join('{jx}*f[{j}][i]'.format(j=j, jx=jx)
                         for j, jx in enumerate(wts) if jx != 0)};

    min_f = f[0][i];
    % for j in range(1, nupts):
    min_f = fmin(min_f, f[${j}][i]);
    % endfor

    if (min_f < ${eps}) {
        beta = abs(avg_f/max(avg_f - min_f, ${eps}));
        beta = fmax(0.0, fmin(beta, 1.0));

        % for j in range(nupts):
        f[${j}][i] = beta*(f[${j}][i] - avg_f) + avg_f;
        % endfor
    }
}

</%pyfr:kernel>
