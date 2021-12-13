# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='calcentropy' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ent_min='out fpdtype_t'>

<% large_number = 10**10 %>

ent_min = ${large_number};
fpdtype_t d,p,e;

% for i in range(nupts):
    d = u[${i}][0];
    % if ndims == 2:
        p = ${c['gamma'] - 1}*(u[${i}][${nvars - 1}] - 
            (0.5/d)*(u[${i}][1]*u[${i}][1] + 
                     u[${i}][2]*u[${i}][2]));
    % elif ndims == 3:
        p = ${c['gamma'] - 1}*(u[${i}][${nvars - 1}] - 
            (0.5/d)*(u[${i}][1]*u[${i}][1] + 
                     u[${i}][2]*u[${i}][2] + 
                     u[${i}][3]*u[${i}][3]));
    % endif

    e = (d <= 0 || p <= 0) ? ${large_number} : d*log(p/pow(d, ${c['gamma']}));

    ent_min = fmin(ent_min, e);
% endfor


</%pyfr:kernel>
