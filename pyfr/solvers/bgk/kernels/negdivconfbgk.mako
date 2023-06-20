# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.dvm'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              f='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              u='in broadcast fpdtype_t[${str(nvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nvars)}]'>
    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};

    for (int j = 0; j < ${nvars}; j++)
    {
        w[0] += f[j]*M[0][j];

        % if ndims == 2:
            w[1] += f[j]*M[0][j]*u[j][0];
            w[2] += f[j]*M[0][j]*u[j][1];

            % if delta:
                w[3] += f[j]*M[0][j]*(0.5*(u[j][0]*u[j][0] + u[j][1]*u[j][1]) + u[j][2]);
            % else:
                w[3] += f[j]*M[0][j]*0.5*(u[j][0]*u[j][0] + u[j][1]*u[j][1]);
            % endif
        % elif ndims == 3:
            w[1] += f[j]*M[0][j]*u[j][0];
            w[2] += f[j]*M[0][j]*u[j][1];
            w[3] += f[j]*M[0][j]*u[j][2];

            % if delta:
                w[4] += f[j]*M[0][j]*(0.5*(u[j][0]*u[j][0] + u[j][1]*u[j][1] + u[j][2]*u[j][2]) + u[j][3]);
            % else:
                w[4] += f[j]*M[0][j]*0.5*(u[j][0]*u[j][0] + u[j][1]*u[j][1] + u[j][2]*u[j][2]);
            % endif
        % endif

    }

    // Convert to primitive
    fpdtype_t rho, U, V, p;
% if ndims == 2:
    rho = w[0];
    U = w[1]/w[0];
    V = w[2]/w[0];
    p = ${c['gamma'] - 1.0}*(w[3] - 0.5*rho*(U*U + V*V));
% elif ndims == 3:
    fpdtype_t W;
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

    ${pyfr.expand('iterate_DVM', 'alpha', 'w', 'u', 'M')};

for (int i = 0; i < ${nvars}; i++)
{
    % if ndims == 2:
        g = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])));
        % if delta:
        fpdtype_t theta = 1.0/(2.0*alpha[1]);
        g *= ${lam}*pow(u[i][2]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][2]/theta);
        % endif
    % elif ndims == 3:
        g = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4])));
        % if delta:
        fpdtype_t theta = 1.0/(2.0*alpha[1]);
        g *= ${lam}*pow(u[i][3]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][3]/theta);
        % endif
    % endif

    tdivtconf[i] = -rcpdjac*tdivtconf[i] + (g - f[i])/${tau};
}


</%pyfr:kernel>
