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

    // Compute right state
    fpdtype_t rhor, Ur[${ndims}], pr;
    fpdtype_t wr[${ndims+2}];
    rhor = rhol;
    
    wr[0] = rhor;
    % for i, v in enumerate('uvw'[:ndims]):
    Ur[${i}] = ${c[v]};
    wr[${i+1}] = rhor*Ur[${i}];
    % endfor
    pr = ${c['theta']}*rhor;
    wr[${ndims+1}] = ${1.0/(c['gamma'] - 1.0)}*pr  + 0.5*rhor*${pyfr.dot('Ur[{i}]', i=ndims)};
    
    // Compute discrete wall Maxwellian
    fpdtype_t Mw[${nvars}];
    
    fpdtype_t alpha[${ndims + 2}];
    alpha[0] = rhor*${(2*pi*c['theta'])**(-ndims/2.0)};
    alpha[1] = ${1.0/(2*c['theta'])};
    alpha[2] = Ur[0];
    alpha[3] = Ur[1];
    % if ndims == 3:
    alpha[4] = Ur[2];
    % endif

    ${pyfr.expand('iterate_DVM', 'alpha', 'wr', 'u', 'M')};

    // Compute mass-preserving scaling factor
    fpdtype_t eta1 = 0.0;
    fpdtype_t eta2 = 0.0;
    fpdtype_t un;

    for (int i = 0; i < ${nvars}; i++)
    {
        un = ${pyfr.dot('u[i][{j}]', 'nl[{j}]', j=ndims)};
        
        % if ndims == 2:
            Mw[i] = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])));
            % if delta:
            fpdtype_t theta = 1.0/(2.0*alpha[1]);
            Mw[i] *= ${lam}*pow(u[i][2]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][2]/theta);
            % endif
        % elif ndims == 3:
            Mw[i] = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4])));
            % if delta:
            fpdtype_t theta = 1.0/(2.0*alpha[1]);
            Mw[i] *= ${lam}*pow(u[i][3]/theta, ${0.5*delta - 1.})*(1./theta)*exp(-u[i][3]/theta);
            % endif
        % endif

        if (un > 0.0) {
            eta1 += fl[i]*M[0][i]*abs(un);
        }
        else {

            eta2 += Mw[i]*M[0][i]*abs(un);
        }
    }

    fpdtype_t eta = eta1/eta2;

    for (int i = 0; i < ${nvars}; i++)
    {
        fr[i] = eta*Mw[i];
    }

</%pyfr:macro>
