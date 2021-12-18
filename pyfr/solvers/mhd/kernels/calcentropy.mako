# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='calcentropy' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin='out fpdtype_t'
              entmin_cpy='out fpdtype_t'>

<% large_number = 10**10 %>

entmin = ${large_number};
fpdtype_t d,p,e;
fpdtype_t rhou,rhov,rhow,E,Bx,By,Bz,BdotB2;

% for i in range(nupts):
    d = u[${i}][0];
    % if ndims == 2:
        rhou = u[${i}][1]; rhov = u[${i}][2];
        Bx = u[${i}][3]; By = u[${i}][4];
        E = u[${i}][6];

        BdotB2 = 0.5*(Bx*Bx + By*By);

        // Compute the pressure
        p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov)
                                 - BdotB2);
    % elif ndims == 3:
        rhou = u[${i}][1]; rhov = u[${i}][2]; rhow = u[${i}][3];
        Bx = u[${i}][4]; By = u[${i}][5]; Bz = u[${i}][6];
        E = u[${i}][8];

        BdotB2 = 0.5*(Bx*Bx + By*By + Bz*Bz);

        // Compute the pressure
        p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov + rhow*rhow)
                                 - BdotB2);
    % endif

    e = (d <= 0 || p <= 0) ? ${large_number} : d*log(p/pow(d, ${c['gamma']}));

    entmin = fmin(entmin, e);
% endfor

entmin_cpy = entmin;
</%pyfr:kernel>
