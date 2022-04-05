# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% large_number = 10**10 %>

<%pyfr:macro name='filter' params='newsolmodes, filtsol, zeta, entmin, withinbounds'>

    dmin = ${large_number}; pmin = ${large_number}; emin = ${large_number};
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        % for i in range(nupts):
            filtmodes[${i}][vidx] = pow(${ffac[i]}, zeta)*newsolmodes[${i}][vidx];
        % endfor
    }

    for (int uidx = 0; uidx < ${nupts}; uidx++) {
        for (int vidx = 0; vidx < ${nvars}; vidx++) {
            filtsol[uidx][vidx] = 0.0;
            for (int midx = 0; midx < ${nupts}; midx++) {
                filtsol[uidx][vidx] += vdm[uidx][midx]*filtmodes[midx][vidx];
            }
        }

        d = filtsol[uidx][0];
        % if ndims == 2:
            p = ${c['gamma'] - 1}*(filtsol[uidx][${nvars - 1}] - 
                (0.5/d)*(filtsol[uidx][1]*filtsol[uidx][1] + 
                         filtsol[uidx][2]*filtsol[uidx][2]));
        % elif ndims == 3:
            p = ${c['gamma'] - 1}*(filtsol[uidx][${nvars - 1}] - 
                (0.5/d)*(filtsol[uidx][1]*filtsol[uidx][1] + 
                         filtsol[uidx][2]*filtsol[uidx][2] + 
                         filtsol[uidx][3]*filtsol[uidx][3]));
        % endif


        e = (d <= ${dtol} || p <= ${ptol}) ? -${large_number} : d*log(p/pow(d, ${c['gamma']}));

        dmin = fmin(dmin, d);
        pmin = fmin(pmin, p);
        emin = fmin(emin, e);

    }

    if (dmin >= ${dtol} && pmin >= ${ptol} && emin >= entmin - ${etol}) {
        withinbounds = 1; 
    }
    else {
        withinbounds = 0;
    }
</%pyfr:macro>


<%pyfr:kernel name='filter' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf_inv='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              tdivtconf_vis='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              entmin='in fpdtype_t'
              vdm='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              invvdm='in fpdtype_t[${str(nupts)}][${str(nupts)}]'>

fpdtype_t newsol[${nupts}][${nvars}];
fpdtype_t du_vis[${nupts}][${nvars}];
fpdtype_t newsolmodes[${nupts}][${nvars}];
fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];

// Compute -divF and forward Euler prediction of next time step
// Separate viscous and inviscid components
for (int uidx = 0; uidx < ${nupts}; uidx++) {
    % for v, ex in enumerate(srcex):
        newsol[uidx][${v}] = u[uidx][${v}] + ${dt}*tdivtconf_inv[uidx][${v}];
        du_vis[uidx][${v}] = ${dt}*(tdivtconf_vis[uidx][${v}] - tdivtconf_inv[uidx][${v}]);
    % endfor
}

// Get modal form at next time step
for (int uidx = 0; uidx < ${nupts}; uidx++) {
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        newsolmodes[uidx][vidx] = 0.0;
        for (int midx = 0; midx < ${nupts}; midx++) {
            newsolmodes[uidx][vidx] += invvdm[uidx][midx]*newsol[midx][vidx];
        }
    }
}



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

for (int uidx = 0; uidx < ${nupts}; uidx++) {
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        newsol[uidx][vidx] = filtsol[uidx][vidx] + du_vis[uidx][vidx];
    }
}

// Get modal form of filtered viscous solution
for (int uidx = 0; uidx < ${nupts}; uidx++) {
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        newsolmodes[uidx][vidx] = 0.0;
        for (int midx = 0; midx < ${nupts}; midx++) {
            newsolmodes[uidx][vidx] += invvdm[uidx][midx]*newsol[midx][vidx];
        }

    }
}


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

for (int uidx = 0; uidx < ${nupts}; uidx++) {
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        tdivtconf_vis[uidx][vidx] = (newsol[uidx][vidx] - u[uidx][vidx])/${dt}; // Store in upts_outb
    }
}

</%pyfr:kernel>


