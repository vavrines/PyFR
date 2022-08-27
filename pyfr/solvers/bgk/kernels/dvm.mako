# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

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
            % if ndims == 2:
                gm = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])));
            % elif ndims == 3:
                gm = alpha[0]*exp(-alpha[1]*((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4])));
            % endif

            // Integrate
            mmnts[0] = 1.0;
            % if ndims == 2:
            mmnts[1] = u[i][0];
            mmnts[2] = u[i][1];
            mmnts[3] = 0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
            % elif ndims == 3:
            mmnts[1] = u[i][0];
            mmnts[2] = u[i][1];
            mmnts[3] = u[i][2];
            mmnts[4] = 0.5*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
            % endif


            % for ivar in range(ndims+2):
                R[${ivar}] += mmnts[${ivar}]*M[0][i]*gm;

                J[${ivar}][0] += mmnts[${ivar}]*M[0][i]*gm/alpha[0];
                J[${ivar}][2] += mmnts[${ivar}]*M[0][i]*gm*2*alpha[1]*(u[i][0] - alpha[2]);
                J[${ivar}][3] += mmnts[${ivar}]*M[0][i]*gm*2*alpha[1]*(u[i][1] - alpha[3]);

                % if ndims == 2:
                    J[${ivar}][1] += mmnts[${ivar}]*M[0][i]*gm*(-((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])));
                % elif ndims == 3:
                    J[${ivar}][1] += mmnts[${ivar}]*M[0][i]*gm*(-((u[i][0]-alpha[2])*(u[i][0]-alpha[2]) + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) + (u[i][2]-alpha[4])*(u[i][2]-alpha[4])));
                    J[${ivar}][4] += mmnts[${ivar}]*M[0][i]*gm*2*alpha[1]*(u[i][2] - alpha[4]);
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