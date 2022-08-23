# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

<%pyfr:macro name='iterate_DVM' params='alpha'>
    fpdtype_t newalpha[${ndims+2}];
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    
    for (int iter = 0; iter < ${niters}; iter++) {
        // Zero cost-function and Jacobian
        % for ivar in range(ndims+2):
            R[${ivar}] = 0; 
            % for jvar in range(ndims+2):
                J[${ivar}][${jvar}] = 0; 
            % endfor
        % endfor

        % for i in range(nvars):
            // Compute discrete Maxwellian
            % if ndims == 2:
                g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0)));
            % elif ndims == 3:
                g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0) + pow(${u[i,2]} - alpha[4], 2.0)));
            % endif

            // Integrate
            % for ivar in range(ndims+2):
                R[${ivar}] += ${PSint[i]*moments[i][ivar]}*g;

                J[${ivar}][0] += ${PSint[i]*moments[i][ivar]}*g/alpha[0];
                J[${ivar}][2] += ${PSint[i]*moments[i][ivar]}*g*2*alpha[1]*(${u[i,0]} - alpha[2]);
                J[${ivar}][3] += ${PSint[i]*moments[i][ivar]}*g*2*alpha[1]*(${u[i,1]} - alpha[3]);

                % if ndims == 2:
                    J[${ivar}][1] += ${PSint[i]*moments[i][ivar]}*g*(-(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0)));
                % elif ndims == 3:
                    J[${ivar}][1] += ${PSint[i]*moments[i][ivar]}*g*(-(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0) + pow(${u[i,2]} - alpha[4], 2.0)));
                    [${ivar}][4] += ${PSint[i]*moments[i][ivar]}*g*2*alpha[1]*(${u[i,2]} - alpha[4]);
                % endif
            % endfor
        % endfor

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
            newalpha[${var}] = alpha[${var}] - (${' + '.join('Jinv[{var}][{i}]*R[{i}]'.format(var=var, i=i) for i in range(ndims+2))});
        % endfor
        
        % for var in range(ndims+2):
            alpha[${var}] = newalpha[${var}];
        % endfor
    }
</%pyfr:macro>