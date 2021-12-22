# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% large_number = 10**10 %>

<%pyfr:macro name='filter' params='newsolmodes, filtsol, zeta, entmin, withinbounds'>

    dmin = ${large_number}; pmin = ${large_number}; emin = ${large_number};
    % for i,v in pyfr.ndrange(nupts, nvars):
        filtmodes[${i}][${v}] = pow(${ffac[i]}, zeta)*newsolmodes[${i}][${v}];
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


<%pyfr:kernel name='negdivconffilter' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf_inv='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              tdivtconf_vis='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              entmin='in fpdtype_t'>

fpdtype_t newsol[${nupts}][${nvars}];
fpdtype_t du_vis[${nupts}][${nvars}];
fpdtype_t newsolmodes[${nupts}][${nvars}];
fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];

// Compute -divF and forward Euler prediction of next time step
// Separate viscous and inviscid components
% for i in range(nupts):
    % for v, ex in enumerate(srcex):
        tdivtconf_inv[${i}][${v}] = -rcpdjac[${i}]*tdivtconf_inv[${i}][${v}] + ${ex};
        tdivtconf_vis[${i}][${v}] = -rcpdjac[${i}]*tdivtconf_vis[${i}][${v}] + ${ex};
        newsol[${i}][${v}] = u[${i}][${v}] + ${dt}*tdivtconf_inv[${i}][${v}];
        du_vis[${i}][${v}] = ${dt}*(tdivtconf_vis[${i}][${v}] - tdivtconf_inv[${i}][${v}]);
    % endfor
% endfor

// Get modal form at next time step
% for i,v in pyfr.ndrange(nupts, nvars):
    newsolmodes[${i}][${v}] = ${' + '.join('{jx}*newsol[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor


// Compute zeta for positivity-preserving and minimum entropy principle satisfying
// ***********************************************************
fpdtype_t zeta_low = 0;
fpdtype_t zeta_high = 0.5;
fpdtype_t zeta = 0;
fpdtype_t pmin, dmin, emin, p, d, e, withinbounds;
fpdtype_t neginf = -${large_number}; 

// Check if solution is already within bounds
${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_low', 'entmin', 'withinbounds')};


// If within bounds, return unfiltered solution
if (withinbounds == 1){
    zeta = 0;
}
// Else apply filter
else {
    // Check that upper bound on filter strength satisfies bounds
    ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_high', 'entmin', 'withinbounds')};

    // If not, increase upper bound
    if (withinbounds == 0){ 
        zeta_high = 12; 
    }

    // Perform bisection method to find zeta
    for (int iter = 0; iter < ${niters}; iter++) {
        zeta = 0.5*(zeta_low + zeta_high);
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta', 'entmin', 'withinbounds')};

        if (withinbounds == 1) {
            zeta_high = zeta; 
        }
        else {
            zeta_low = zeta;
        }
    }

    // Take bounds-preserving value for zeta and compute filtered solution
    zeta = zeta_high;
    ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta', 'entmin', 'withinbounds')};
}


// ***********************************************************
// Add viscous component and check if positivity-preserving
% for i,v in pyfr.ndrange(nupts, nvars):
    newsol[${i}][${v}] = filtsol[${i}][${v}] + du_vis[${i}][${v}];
% endfor

// Get modal form of filtered viscous solution
% for i,v in pyfr.ndrange(nupts, nvars):
    newsolmodes[${i}][${v}] = ${' + '.join('{jx}*newsol[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor

${pyfr.expand('filter', 'newsolmodes', 'newsol', 'zeta', 'neginf', 'withinbounds')};

// If not positivity-preserving, continue increasing filter strength
if (withinbounds == 0) {
    zeta_low = zeta;
    zeta_high = 12;

    for (int iter = 0; iter < ${niters}; iter++) {
        zeta = 0.5*(zeta_low + zeta_high);
        ${pyfr.expand('filter', 'newsolmodes', 'newsol', 'zeta', 'neginf', 'withinbounds')};

        if (withinbounds == 1) {
            zeta_high = zeta; 
        }
        else {
            zeta_low = zeta;
        }
    }
    zeta = zeta_high;
    
    ${pyfr.expand('filter', 'newsolmodes', 'newsol', 'zeta', 'neginf', 'withinbounds')};
}


// Compute forward Euler approximation of -divF
% for i,v in pyfr.ndrange(nupts, nvars):
    tdivtconf_inv[${i}][${v}] = (newsol[${i}][${v}] - u[${i}][${v}])/${dt}; // Store in upts_outb
% endfor

</%pyfr:kernel>


