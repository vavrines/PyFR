# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.bgk.kernels.dvm'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    // Compute left conservative state
    fpdtype_t wl[${ndims+2}] = {0};
    for (int i = 0; i < ${nvars}; i++)
    {
        wl[0] += fl[i]*M[0][i];
        wl[1] += fl[i]*M[0][i]*u[i][0];
        wl[2] += fl[i]*M[0][i]*u[i][1];

        % if ndims == 2:
            % if delta:
                wl[3] += fl[i]*M[0][i]*(0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1]) + u[i][2]);
            % else:
                wl[3] += fl[i]*M[0][i]*0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
            % endif
        % elif ndims == 3:
            wl[3] += fl[i]*M[0][i]*u[i][2];            
            % if delta:
                wl[4] += fl[i]*M[0][i]*(0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]) + u[i][3]);
            % else:
                wl[4] += fl[i]*M[0][i]*0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
            % endif
        % endif
    }

    // Compute primitives
    fpdtype_t rhol, Ul[${ndims}], pl;

    rhol = wl[0];
    fpdtype_t invrhol = 1.0/rhol;
    % for i in range(ndims):
    Ul[${i}] = invrhol*wl[${i + 1}];
    % endfor
    pl = ${c['gamma'] - 1.0}*(wl[${ndims+1}] - 0.5*rhol*${pyfr.dot('Ul[{i}]', i=ndims)});

    // Compute RHS state
    fpdtype_t w[${ndims + 2}];
    w[0] = rhol;
% for i in range(ndims):
    w[${i + 1}] = wl[${i + 1}];
% endfor
    w[${ndims+1}] = ${c['p']}/${(c['gamma'] - 1)}
                    + 0.5*(1.0/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims + 1))};
    
    fpdtype_t rho = rhol;
    fpdtype_t p = ${c['p']};
    fpdtype_t theta = p/rho;
    
    fpdtype_t alpha[${ndims + 2}];
    alpha[0] = rho*pow(${2*pi}*theta, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta);
    alpha[2] = Ul[0];
    alpha[3] = Ul[1];
    % if ndims == 3:
    alpha[4] = Ul[2];
    % endif

    ${pyfr.expand('iterate_DVM', 'alpha', 'w', 'u', 'M')};

    for (int i = 0; i < ${nvars}; i++)
    {
        % if ndims == 2:
            fr[i] = alpha[0]*exp(-alpha[1]*(pow(u[i][0] - alpha[2], 2.0) + pow(u[i][1] - alpha[3], 2.0)));
            % if delta:
            fpdtype_t theta = 1.0/(2.0*alpha[1]);
            fr[i] *= ${lam}*pow(u[i][2]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][2]/theta);
            % endif
        % elif ndims == 3:
            fr[i] = alpha[0]*exp(-alpha[1]*(pow(u[i][0] - alpha[2], 2.0) + pow(u[i][1] - alpha[3], 2.0) + pow(u[i][2] - alpha[4], 2.0)));
            % if delta:
            fpdtype_t theta = 1.0/(2.0*alpha[1]);
            fr[i] *= ${lam}*pow(u[i][3]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][3]/theta);
            % endif
        % endif
    }
    
</%pyfr:macro>
