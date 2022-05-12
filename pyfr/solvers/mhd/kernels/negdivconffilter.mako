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
            rhou = filtsol[uidx][1]; rhov = filtsol[uidx][2];
            Bx = filtsol[uidx][3]; By = filtsol[uidx][4];
            E = filtsol[uidx][6];

            BdotB2 = 0.5*(Bx*Bx + By*By);

            // Compute the pressure
            p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov)
                                     - BdotB2);
        % elif ndims == 3:
            rhou = filtsol[uidx][1]; rhov = filtsol[uidx][2]; rhow = filtsol[uidx][3];
            Bx = filtsol[uidx][4]; By = filtsol[uidx][5]; Bz = filtsol[uidx][6];
            E = filtsol[uidx][8];

            BdotB2 = 0.5*(Bx*Bx + By*By + Bz*Bz);

            // Compute the pressure
            p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov + rhow*rhow)
                                     - BdotB2);
        % endif


        e = (d <= ${dtol} || p <= ${ptol}) ? ${large_number} : d*log(p/pow(d, ${c['gamma']}));

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


<%pyfr:kernel name='negdivconffilter' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              entmin='in fpdtype_t'
              vdm='in fpdtype_t[${str(nupts)}][${str(nupts)}]'
              invvdm='in fpdtype_t[${str(nupts)}][${str(nupts)}]'>

fpdtype_t newsol[${nupts}][${nvars}];
fpdtype_t newsolmodes[${nupts}][${nvars}];
fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];
fpdtype_t source[${nupts}][${nvars}];


fpdtype_t rcprho, v[${ndims}], B[${ndims}], Bdotv, divB;
fpdtype_t rhou, rhov, rhow, E, Bx, By, Bz, BdotB2;

// Compute -divF and forward Euler prediction of next time step
fpdtype_t divf;
for (int uidx = 0; uidx < ${nupts}; uidx++) {
    rcprho = 1.0/u[uidx][0];

    % if ndims == 2:
        // Velocity and magnetic fields
        v[0] = rcprho*u[uidx][1]; v[1] = rcprho*u[uidx][2];
        B[0] = u[uidx][3]; B[1] = u[uidx][4];

        // Compute B·v
        Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=2)};
        divB = tdivtconf[uidx][5];
        u[uidx][5] = divB;
        tdivtconf[uidx][5] = 0.0;

        // Untransform the divergences and apply the source terms
        // Density
        source[uidx][0] = 0.0;
        // Momentum
        source[uidx][1] = -rcpdjac[uidx]*divB*B[0] + ${srcex[1]};
        source[uidx][2] = -rcpdjac[uidx]*divB*B[1] + ${srcex[2]};
        // Magnetic field
        source[uidx][3] = -rcpdjac[uidx]*divB*v[0] + ${srcex[3]};
        source[uidx][4] = -rcpdjac[uidx]*divB*v[1] + ${srcex[4]};
        // DivB
        source[uidx][5] = 0.0; 
        // Energy
        source[uidx][6] = -rcpdjac[uidx]*divB*Bdotv + ${srcex[6]}; 

    % elif ndims == 3:
        // Velocity and magnetic fields
        v[0] = rcprho*u[uidx][1]; v[1] = rcprho*u[uidx][2]; v[2] = rcprho*u[uidx][3];
        B[0] = u[uidx][4]; B[1] = u[uidx][5]; B[2] = u[uidx][6];

        // Compute B·v
        Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=3)};
        divB = tdivtconf[uidx][7];
        u[uidx][7] = divB;
        tdivtconf[uidx][7] = 0.0; 

        // Untransform the divergences and apply the source terms
        // Density
        source[uidx][0] = 0.0 + ${srcex[0]};
        // Momentum
        source[uidx][1] = -rcpdjac[uidx]*divB*B[0] + ${srcex[1]};
        source[uidx][2] = -rcpdjac[uidx]*divB*B[1] + ${srcex[2]};
        source[uidx][3] = -rcpdjac[uidx]*divB*B[2] + ${srcex[3]};
        // Magnetic field
        source[uidx][4] = -rcpdjac[uidx]*divB*v[0] + ${srcex[4]};
        source[uidx][5] = -rcpdjac[uidx]*divB*v[1] + ${srcex[5]};
        source[uidx][6] = -rcpdjac[uidx]*divB*v[2] + ${srcex[6]};
        // DivB
        source[uidx][7] = 0.0; 
        // Energy
        source[uidx][8] = -rcpdjac[uidx]*divB*Bdotv + ${srcex[8]}; 
    % endif

    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        tdivtconf[uidx][vidx] = -rcpdjac[uidx]*tdivtconf[uidx][vidx];
    }

    % for v, ex in enumerate(srcex):
        newsol[uidx][${v}] = u[uidx][${v}] + ${dt}*tdivtconf[uidx][${v}];
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



for (int uidx = 0; uidx < ${nupts}; uidx++) {
    for (int vidx = 0; vidx < ${nvars}; vidx++) {
        tdivtconf[uidx][vidx] = (filtsol[uidx][vidx] - u[uidx][vidx])/${dt} + source[uidx][vidx]; // Store in upts_outb
    }
    
    % if ndims == 2:
        tdivtconf[uidx][5] = 0.0;
    % elif ndims == 3:
        tdivtconf[uidx][7] = 0.0; 
    % endif
}

</%pyfr:kernel>
