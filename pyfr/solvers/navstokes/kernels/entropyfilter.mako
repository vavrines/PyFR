# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% large_number = 10**10 %>

<%pyfr:macro name='filter_ent' params='u_modes, filtsol, zeta, entmin, withinbounds'>

    dmin = ${large_number}; pmin = ${large_number}; emin = ${large_number};
    % for i,v in pyfr.ndrange(nupts, nvars):
        filtmodes[${i}][${v}] = pow(${ffac[i]}, zeta)*u_modes[${i}][${v}];
    % endfor

    % for i in range(nupts):
        % for v in range(nvars):
            filtsol[${i}][${v}] = ${' + '.join('{jx}*filtmodes[{j}][{v}]'.format(j=j, jx=jx, v=v) 
                                                for j, jx in enumerate(vdm[i]) if jx != 0)};
        % endfor

        d = filtsol[${i}][0];
        % if ndims == 2:
            p = ${c['gamma'] - 1}*(filtsol[${i}][${nvars - 1}] - 
                (0.5/d)*(filtsol[${i}][1]*filtsol[${i}][1] + 
                         filtsol[${i}][2]*filtsol[${i}][2]));
        % elif ndims == 3:
            p = ${c['gamma'] - 1}*(filtsol[${i}][${nvars - 1}] - 
                (0.5/d)*(filtsol[${i}][1]*filtsol[${i}][1] + 
                         filtsol[${i}][2]*filtsol[${i}][2] + 
                         filtsol[${i}][3]*filtsol[${i}][3]));
        % endif

        e = (d <= ${dtol} || p <= ${ptol}) ? -${large_number} : d*log(p/pow(d, ${c['gamma']}));

        dmin = fmin(dmin, d);
        pmin = fmin(pmin, p);
        emin = fmin(emin, e);
    % endfor

    if (dmin >= ${dtol} && pmin >= ${ptol} && emin >= entmin - ${etol}) {
        withinbounds = 1; 
    }
    else {
        withinbounds = 0;
    }
</%pyfr:macro>


<%pyfr:kernel name='entropyfilter' ndim='1'
              t='scalar fpdtype_t'
              u_modes='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin='in fpdtype_t'>

fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];

fpdtype_t zeta_low = 0;
fpdtype_t zeta_high = 0.5;
fpdtype_t zeta = 0;
fpdtype_t pmin, dmin, emin, p, d, e, withinbounds;

// Compute zeta for positivity-preserving and minimum entropy principle satisfying

// Check if solution is already within bounds
${pyfr.expand('filter_ent', 'u_modes', 'filtsol', 'zeta_low', 'entmin', 'withinbounds')};

// If within bounds, return unfiltered solution, else apply filter
if (withinbounds == 0){
    // Check that upper bound on filter strength satisfies bounds
    ${pyfr.expand('filter_ent', 'u_modes', 'filtsol', 'zeta_high', 'entmin', 'withinbounds')};

    // If not, increase upper bound
    if (withinbounds == 0){ 
        zeta_high = 12; 
    }

    // Perform bisection method to find zeta
    for (int iter = 0; iter < ${niters}; iter++) {
        zeta = 0.5*(zeta_low + zeta_high);
        ${pyfr.expand('filter_ent', 'u_modes', 'filtsol', 'zeta', 'entmin', 'withinbounds')};

        if (withinbounds == 1) {
            zeta_high = zeta; 
        }
        else {
            zeta_low = zeta;
        }
    }

    // Take bounds-preserving value for zeta and compute filtered solution
    zeta = zeta_high;

    ${pyfr.expand('filter_ent', 'u_modes', 'filtsol', 'zeta', 'entmin', 'withinbounds')};

}

// Store filtered inviscid solution in place of modal form
% for i,v in pyfr.ndrange(nupts, nvars):
    u_modes[${i}][${v}] = filtsol[${i}][${v}];
% endfor

</%pyfr:kernel>


