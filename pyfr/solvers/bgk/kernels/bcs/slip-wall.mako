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
            wl[3] += fl[i]*M[0][i]*0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
        % elif ndims == 3:
            wl[3] += fl[i]*M[0][i]*u[i][2];
            lw[4] += fl[i]*M[0][i]*0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
        % endif
    }

    // Compute primitives
    fpdtype_t rhol, Ul[${ndims}], pl;

    rhol = wl[0];
    fpdtype_t invrhol = 1.0/rhol;
    % for i in range(ndims):
    Ul[${i}] = invrhol*wl[${i + 1}];
    % endfor
    pl = ${c['gamma'] - 1}*(wl[${ndims+1}] - 0.5*rhol*${pyfr.dot('Ul[{i}]', i=ndims)});

    // Compute right state (subtract normal velocity)
    fpdtype_t rhor, Ur[${ndims}], pr;
    fpdtype_t wr[${ndims+2}];
    wr[0] = rhor = rhol;
    fpdtype_t nor = ${' + '.join(f'Ul[{i}]*nl[{i}]' for i in range(ndims))};
    % for i in range(ndims):
    Ur[${i}] = Ul[${i}] - nor*nl[${i}];
    wr[${i+1}] = rhor*Ur[${i}];
    % endfor
    pr = pl;
    wr[${ndims+1}] = pr/(${c['gamma'] - 1}) + 0.5*rhor*${pyfr.dot('Ur[{i}]', i=ndims)};
    
    // Compute discrete wall Maxwellian
    fpdtype_t Mw[${nvars}];
    fpdtype_t theta = pr/rhor;
    
    fpdtype_t alpha[${ndims + 2}];
    alpha[0] = rhor*pow(${2*pi}*theta, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta);
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
        % elif ndims == 3:
        Mw[i] = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4])));
        % endif

        if (un > 0.0) {
            eta1 += fl[i]*M[0][i];
        }
        else {

            eta2 += Mw[i]*M[0][i];
        }
    }

    fpdtype_t eta = eta1/max(eta2, ${10**-12});

    for (int i = 0; i < ${nvars}; i++)
    {
        fr[i] = eta*Mw[i];
    }

</%pyfr:macro>
