# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='normalizeresidual' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              r='inout fpdtype_t[${str(nupts)}]'
              usmats='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              >


// Calculate entropy
% for i in range(nupts):
    u[${i}][${nvars-1}] = 0.0;
% endfor

// Calculate rho, P, entropy
fpdtype_t ub[${nupts}][${kvars}];
fpdtype_t dsdr[${nupts}], dsdp[${nupts}];
fpdtype_t tmp;
% for i in range(nupts):
    ub[${i}][0] = u[${i}][0]; // Rho
    ub[${i}][1] = ${c['gamma'] - 1}*(u[${i}][${nvars-2}] - 0.5*(pow(u[${i}][1], 2.0) + pow(u[${i}][2], 2.0))/u[${i}][0]); // Pressure

    // Entropy
    u[${i}][${nvars-1}] = -ub[${i}][0]*log(ub[${i}][1]*pow(ub[${i}][0], -${c['gamma']}))/(${c['gamma']- 1.0});
    ub[${i}][2] = u[${i}][${nvars-1}];

    dsdr[${i}] = -ub[${i}][2]/(${c['gamma'] - 1}) + ${c['gamma']/(c['gamma'] - 1)};
    dsdp[${i}] = -u[${i}][0]/(${c['gamma'] - 1}*ub[${i}][1]);
% endfor

fpdtype_t gradu[${nupts}][${ndims}][${kvars}];
fpdtype_t tgradu[${nupts}][${ndims}][${kvars}];

fpdtype_t line_sol[${order+1}][${kvars}];

// Perform along xi direction
% for j in range(order+1):
    // Gather solution along constant eta line
    % for i in range(order+1):
        % for var in range(kvars):
            line_sol[${i}][${var}] = ub[${i+ j*(order+1)}][${var}];
        % endfor
    % endfor

    // Calculate du/dxi at solution points
    % for var in range(kvars):
        % for i in range(order+1):
            tmp = ${' + '.join('{mx}*line_sol[{m}][{svar}]'.format(m=m, mx=mx, svar=var)
                       for m, mx in enumerate(diffmatSOL[i]) if mx != 0)};
            gradu[${i+ j*(order+1)}][0][${var}] = tmp;
        % endfor
    % endfor 
% endfor

// Perform along eta direction
% for i in range(order+1):
    // Gather solution along constant xi line
    % for j in range(order+1):
        % for var in range(kvars):
            line_sol[${j}][${var}] = ub[${i+ j*(order+1)}][${var}];
        % endfor
    % endfor

    // Calculate du/deta at solution points
    % for var in range(kvars):
        % for j in range(order+1):
            tmp =  ${' + '.join('{mx}*line_sol[{m}][{svar}]'.format(m=m, mx=mx, svar=var)
                       for m, mx in enumerate(diffmatSOL[j]) if mx != 0)};
            gradu[${i+ j*(order+1)}][1][${var}] = tmp;
        % endfor
    % endfor 
% endfor

% for i, dim, var in pyfr.ndrange(nupts, ndims, kvars):
    tgradu[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*gradu[{0}][{2}][{3}]'
                                             .format(i, ndims*dim+k, k, var) for k in range(ndims))};

    tgradu[${i}][${dim}][${var}] = tgradu[${i}][${dim}][${var}]*rcpdjac[${i}];
% endfor

% for i, dim, var in pyfr.ndrange(nupts, ndims, kvars):
    tgradu[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*gradu[{0}][{2}][{3}]'
                                             .format(i, ndims*dim+k, k, var) for k in range(ndims))};

    tgradu[${i}][${dim}][${var}] = tgradu[${i}][${dim}][${var}]*rcpdjac[${i}];
% endfor


fpdtype_t r0 = 0, r1, r2, divs, divr, divp, lambda;
% for i in range(nupts):
    r0 = 0; r1 = 0; r2 = 0; divs = 0; divr = 0; divp = 0;
    % for dim in range(ndims):
        divr += tgradu[${i}][${dim}][0];
        divp += tgradu[${i}][${dim}][1];
        divs += tgradu[${i}][${dim}][2];
    % endfor

    lambda = abs(pow(pow(u[${i}][1]/u[${i}][0], 2.0) + pow(u[${i}][2]/u[${i}][0], 2.0), 0.5)) + sqrt(ub[${i}][1]/ub[${i}][0]);

    r1 = lambda*(abs(divs) + ${tol});
    r2 = lambda*(abs(divr)*abs(dsdr[${i}]) + abs(divp)*abs(dsdp[${i}]) + ${tol});
    r0 = max(r0, max(r1, r2));
% endfor

% for i in range(nupts):
    //r[${i}] = abs(r[${i}])/r0;
% endfor

</%pyfr:kernel>
