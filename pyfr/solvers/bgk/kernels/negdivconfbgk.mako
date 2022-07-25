# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

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

    // Create discrete velocity preserving Maxwellian
    fpdtype_t alpha[${ndims+2}], newalpha[${ndims+2}];
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    fpdtype_t det, g;

    fpdtype_t lam = 0.5*(rho/e);
    alpha[0] = rho*pow(lam/${pi}, ${ndims/2.0});
    alpha[1] = lam;
    alpha[2] = U;
    alpha[3] = V;
% if ndims == 3:
    alpha[4] = W;
% endif 

for (int iter = 0; iter < ${niters}; iter++) {
    
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
        J[${ivar}][4] += ${PSint[i]*moments[i][ivar]}*g*2*alpha[1]*(${u[i,2]} - alpha[4]);
        
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


% for i in range(nvars):
    % if ndims == 2:
    g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0)));
    % elif ndims == 3:
    g = alpha[0]*exp(-alpha[1]*(pow(${u[i,0]} - alpha[2], 2.0) + pow(${u[i,1]} - alpha[3], 2.0) + pow(${u[i,2]} - alpha[4], 2.0)));
    % endif

    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + (g - f[${i}])/${tau};
% endfor

</%pyfr:kernel>
