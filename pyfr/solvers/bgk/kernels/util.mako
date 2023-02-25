# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

<%pyfr:macro name='compute_moments' params='f, u, M, w'>
    fpdtype_t fm;
    for (int i = 0; i < ${nvars}; i++) {
        fm = f[i]*M[0][i];

        w[0] += fm;
        w[1] += fm*u[i][0];
        w[2] += fm*u[i][1];
    % if ndims == 2:
        w[3] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
    % elif ndims == 3:
        w[3] += fm*u[i][2];
        w[4] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
    % endif

    % if delta:
        w[${ndims+1}] += fm*u[i][${ndims}];
    % endif
    }
</%pyfr:macro>

<%pyfr:macro name='con_to_pri' params='w, q'>
    q[0] = w[0];
    q[1] = w[1]/w[0];
    q[2] = w[2]/w[0];
    % if ndims == 2:
    q[3] = ${c['gamma']-1.0}*(w[3] - 0.5*q[0]*(q[1]*q[1] + q[2]*q[2]));
    % elif ndims == 3:
    q[3] = w[3]/w[0];
    q[4] = ${c['gamma']-1.0}*(w[4] - 0.5*q[0]*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]));
    % endif
</%pyfr:macro>

<%pyfr:macro name='pri_to_con' params='q, w'>
    w[0] = q[0];
    w[1] = q[0]*q[1];
    w[2] = q[0]*q[2];
    % if ndims == 2:
    w[3] = ${1.0/(c['gamma']-1.0)}*q[3] + 0.5*q[0]*(q[1]*q[1] + q[2]*q[2]);
    % elif ndims == 3:
    w[3] = q[0]*q[3];
    w[4] = ${1.0/(c['gamma']-1.0)}*q[4] + 0.5*q[0]*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    % endif
</%pyfr:macro>

<%pyfr:macro name='compute_alpha' params='q, alpha'>
    fpdtype_t theta_tmp = q[${ndims+1}]/q[0];
    alpha[0] = q[0]*pow(${2*pi}*theta_tmp, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta_tmp);
    % for i in range(ndims):
    alpha[${i+2}] = q[${i+1}];
    % endfor
</%pyfr:macro>

<%pyfr:macro name='compute_equilibrium_distribution' params='alpha, u, i, g'>
    // Compute square of pecular velocity
    fpdtype_t dv2;
    % if ndims == 2:
    dv2 = (u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]);
    % elif ndims == 3:
    dv2 = (u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4]);
    % endif

    // Compute monatomic Maxwellian
    g = alpha[0]*exp(-alpha[1]*dv2);

    // Amend with internal energy terms if needed
    % if delta:
    // Using theta = 1.0/(2.0*alpha[1])
    g *= ${lam}*pow(2*u[i][${ndims}]*alpha[1], ${0.5*delta - 1.})*(1./theta)*exp(-2*u[i][${ndims}]*alpha[1]);
    % endif
</%pyfr:macro>

<%pyfr:macro name='iterate_DVM' params='alpha, w, u, M'>
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    fpdtype_t mmnts[${ndims+2}];
    fpdtype_t gm;
    
    for (int iter = 0; iter < ${niters}; iter++) {
        // Zero cost-function and Jacobian
        % for ivar in range(ndims+2):
        R[${ivar}] = 0; 
        % for jvar in range(ndims+2):
        J[${ivar}][${jvar}] = 0; 
        % endfor
        % endfor

        // Compute discrete Maxwellian
        for (int i = 0; i < ${nvars}; i++) {
            ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'i', 'gm')};

            // Precompute moment factors
            mmnts[0] = M[0][i]*gm;
            mmnts[1] = M[0][i]*gm*u[i][0];
            mmnts[2] = M[0][i]*gm*u[i][1];
        % if ndims == 2:
            mmnts[3] = 0.5*M[0][i]*gm*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
        % elif ndims == 3:
            mmnts[3] = M[0][i]*gm*u[i][2];
            mmnts[4] = 0.5*M[0][i]*gm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
        % endif

        % if delta:
            mmnts[${ndims+1}] += M[0][i]*gm*u[i][${ndims}];
        % endif 

        % for ivar in range(ndims+2):
            R[${ivar}] += mmnts[${ivar}];

            J[${ivar}][0] += mmnts[${ivar}]/alpha[0];
            J[${ivar}][2] += mmnts[${ivar}]*2*alpha[1]*(u[i][0] - alpha[2]);
            J[${ivar}][3] += mmnts[${ivar}]*2*alpha[1]*(u[i][1] - alpha[3]);

            % if ndims == 2:
            J[${ivar}][1] += -mmnts[${ivar}]*( (u[i][0]-alpha[2])*(u[i][0]-alpha[2])
                                             + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) );
            % elif ndims == 3:
            J[${ivar}][1] += -mmnts[${ivar}]*( (u[i][0]-alpha[2])*(u[i][0]-alpha[2])
                                             + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])
                                             + (u[i][2]-alpha[4])*(u[i][2]-alpha[4]) );
            J[${ivar}][4] += mmnts[${ivar}]*2*alpha[1]*(u[i][2] - alpha[4]);
            % endif

            % if delta:
            J[${ivar}][1] += mmnts[${ivar}]*(${delta} - 4*u[i][${ndims}]*alpha[1])/(2*alpha[1]);
            % endif
        % endfor
        }


        // Get defect
        % for var in range(ndims+2):
        R[${var}] -= w[${var}]; 
        % endfor

        // Compute inverse Jacobian
        % if ndims == 2:
        ${pyfr.expand('compute_4x4inverse', 'J', 'Jinv')};
        % elif ndims == 3:
        ${pyfr.expand('compute_5x5inverse', 'J', 'Jinv')};
        % endif

        // Take Newton iteration
        % for var in range(ndims+2):
        alpha[${var}] = alpha[${var}] - (${' + '.join('Jinv[{var}][{i}]*R[{i}]'.format(var=var, i=i) for i in range(ndims+2))});
        % endfor
    }
</%pyfr:macro>



