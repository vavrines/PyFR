# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='normalizeresidual' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              r='inout fpdtype_t[${str(nupts)}]'
              divu='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
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
% for i in range(nupts):
    ub[${i}][0] = u[${i}][0]; // Rho
    ub[${i}][1] = ${c['gamma'] - 1}*(u[${i}][${nvars-2}] - 0.5*(pow(u[${i}][1], 2.0) + pow(u[${i}][2], 2.0))/u[${i}][0]); // Pressure

    // Entropy
    u[${i}][${nvars-1}] = -ub[${i}][0]*log(ub[${i}][1]*pow(ub[${i}][0], -${c['gamma']}))/(${c['gamma']- 1.0});
    ub[${i}][2] = u[${i}][${nvars-1}];

    dsdr[${i}] = -ub[${i}][2]/(${c['gamma'] - 1}) + ${c['gamma']/(c['gamma'] - 1)};
    dsdp[${i}] = -u[${i}][0]/(${c['gamma'] - 1}*ub[${i}][1]);
% endfor

// divr = 0
// divp = nvars-2
// divs = nvars-1

fpdtype_t r0[${nupts}], smooth_r[${nupts}], smooth_r0[${nupts}], max_r0 = ${tol}, mean_r0 = ${tol}, r1, r2, divs, divr, divp, lambda;
% for i in range(nupts):
    divr = divu[${i}][0];
    divp = divu[${i}][${nvars-2}];
    divs = divu[${i}][${nvars-1}];

    lambda = abs(pow(pow(u[${i}][1]/u[${i}][0], 2.0) + pow(u[${i}][2]/u[${i}][0], 2.0), 0.5)) + sqrt(ub[${i}][1]/ub[${i}][0]);

    r1 = rcpdjac[${i}]*lambda*(abs(divs) + ${tol});
    r2 = rcpdjac[${i}]*lambda*(abs(divr)*abs(dsdr[${i}]) + abs(divp)*abs(dsdp[${i}]) + ${tol});

    r0[${i}] = fmax(${tol}, fmax(r1, r2));
    max_r0 = fmax(max_r0, r0[${i}]);
    mean_r0 = ${quadwts[i]/4.0}*r0[${i}];

    r[${i}] = abs(r[${i}]);
% endfor


% for i in range(nupts):
    smooth_r[${i}] = ${' + '.join('{mx}*r[{m}]'.format(m=m, mx=mx)
                       for m, mx in enumerate(smoothmat[i]) if mx != 0)};
    smooth_r0[${i}] = ${' + '.join('{mx}*r0[{m}]'.format(m=m, mx=mx)
                       for m, mx in enumerate(smoothmat[i]) if mx != 0)};
% endfor


% for i in range(nupts):
    r[${i}] = abs(smooth_r[${i}])/smooth_r0[${i}];
    //r[${i}] = abs(r[${i}]/r0[${i}]);
    //r[${i}] = abs(smooth_r[${i}])/r0[${i}];
    //r[${i}] = abs(smooth_r[${i}])/max_r0;
    //r[${i}] = abs(smooth_r[${i}])/mean_r0;
% endfor

</%pyfr:kernel>
