# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.bgk.kernels.dvm'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    fpdtype_t w[${ndims + 2}];
    w[0] = ${c['rho']};
% for i, v in enumerate('uvw'[:ndims]):
    w[${i + 1}] = (${c['rho']})*(${c[v]});
% endfor
    w[${ndims+1}] = ${c['p']}/${c['gamma'] - 1} +
                       0.5*(1.0/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims + 1))};

    
    fpdtype_t rho = ${c['rho']};
    fpdtype_t p = ${c['p']};
    fpdtype_t theta = p/rho;
    
    fpdtype_t alpha[${ndims + 2}];
    alpha[0] = rho*pow(${2*pi}*theta, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta);
    alpha[2] = ${c['u']};
    alpha[3] = ${c['v']};
    % if ndims == 3:
    alpha[4] = ${c['w']};
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
