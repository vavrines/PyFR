# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.dvm'/>

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
    fpdtype_t rho, U, V, p;
    rho = w[0];
    U = w[1]/w[0];
    V = w[2]/w[0];
    p = ${c['gamma'] - 1.0}*(w[3] - 0.5*rho*(U*U + V*V));
% elif ndims == 3:
    fpdtype_t rho, U, V, W, e;
    rho = w[0];
    U = w[1]/w[0];
    V = w[2]/w[0];
    W = w[3]/w[0];
    p = ${c['gamma'] - 1.0}*(w[4] - 0.5*rho*(U*U + V*V + W*W));
% endif

    // Create discrete velocity preserving Maxwellian
    fpdtype_t alpha[${ndims+2}], g;
    fpdtype_t theta = p/rho;
    alpha[0] = rho*pow(${2*pi}*theta, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta);
    alpha[2] = U;
    alpha[3] = V;
% if ndims == 3:
    alpha[4] = W;
% endif 

    ${pyfr.expand('iterate_DVM', 'alpha')};

% for i in range(nvars):
    % if ndims == 2:
    g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0)));
    % elif ndims == 3:
    g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0) + pow(${u[i,2]} - alpha[4], 2.0)));
    % endif

    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + (g - f[${i}])/${tau};
% endfor

</%pyfr:kernel>
