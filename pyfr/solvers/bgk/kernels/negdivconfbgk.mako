# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              f='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
% for i in range(ndims+2):
    w[${i}] = ${' + '.join('{jx}*f[{j}]*{m}'.format(j=j, jx=jx, m=moments[j][i])
                       for j, jx in enumerate(PSint) if jx != 0)};
% endfor

    // Convert to primitive
% if ndims == 2:
    fpdtype_t rho, U, V, e;
    rho = w[0];
    U = w[1]/w[0];
    V = w[2]/w[0];
    e = w[3] - 0.5*rho*(U*U + V*V);
% elif ndims == 3:
    fpdtype_t rho, U, V, W, e;
    rho = w[0];
    U = w[1]/w[0];
    V = w[2]/w[0];
    W = w[3]/w[0];
    e = w[4] - 0.5*rho*(U*U + V*V + W*W);
% endif

    // Create Maxwellian
    fpdtype_t g, dv2;

% if quasi1d:
    fpdtype_t lam = 0.25*(rho/e);
% else:
    fpdtype_t lam = 0.5*(rho/e);
% endif
% for i in range(nvars):
    % if ndims == 2:
    dv2 = pow(${u[i,0]} - U, 2.0) + pow(${u[i,1]} - V, 2.0);
    % elif ndims == 3:
    dv2 = pow(${u[i,0]} - U, 2.0) + pow(${u[i,1]} - V, 2.0) + pow(${u[i,2]} - W, 2.0);
    % endif

    g = rho*pow(lam/${pi}, ${(quasi1d or ndims)/2.0})*exp(-lam*dv2);

    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + (g - f[${i}])/${tau};
% endfor


</%pyfr:kernel>
