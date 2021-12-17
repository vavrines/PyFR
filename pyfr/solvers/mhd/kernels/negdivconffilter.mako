# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% large_number = 10**10 %>

<%pyfr:macro name='filter' params='newsolmodes, filtsol, zeta, ent_min, withinbounds'>

    dmin = ${large_number}; pmin = ${large_number}; emin = ${large_number};
    % for i,v in pyfr.ndrange(nupts, nvars):
        filtmodes[${i}][${v}] = pow(${ffac[i]}, zeta)*newsolmodes[${i}][${v}];
    % endfor

    % for i,v in pyfr.ndrange(nupts, nvars):
        filtsol[${i}][${v}] = ${' + '.join('{jx}*filtmodes[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(vdm[i]) if jx != 0)};

        d = filtsol[${i}][0];
        % if ndims == 2:
            rhou = u[${i}][1]; rhov = u[${i}][2];
            Bx = u[${i}][3]; By = u[${i}][4];
            E = u[${i}][6];

            BdotB2 = 0.5*(Bx*Bx + By*By);

            // Compute the pressure
            p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov)
                                     - BdotB2);
        % elif ndims == 3:
            rhou = u[${i}][1]; rhov = u[${i}][2]; rhow = u[${i}][3];
            Bx = u[${i}][4]; By = u[${i}][5]; Bz = u[${i}][6];
            E = u[${i}][8];

            BdotB2 = 0.5*(Bx*Bx + By*By + Bz*Bz);

            // Compute the pressure
            p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov + rhow*rhow)
                                     - BdotB2);
        % endif

        e = (d <= ${dtol} || p <= ${ptol}) ? -${large_number} : d*log(p/pow(d, ${c['gamma']}));

        dmin = fmin(dmin, d);
        pmin = fmin(pmin, p);
        emin = fmin(emin, e);
    % endfor

    if (dmin >= ${dtol} && pmin >= ${ptol} && emin >= ent_min - ${etol}) {
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
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              ent_min='in fpdtype_t'>

fpdtype_t newsol[${nupts}][${nvars}];
fpdtype_t newsolmodes[${nupts}][${nvars}];
fpdtype_t filtsol_pp[${nupts}][${nvars}];
fpdtype_t filtsol_mep[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];
fpdtype_t rcprho, v[${ndims}], B[${ndims}], Bdotv, divB;
fpdtype_t rhou, rhov, rhow, E, Bx, By, Bz, BdotB2;

// Compute -divF
% for i in range(nupts):
    rcprho = 1.0/u[${i}][0];
    % if ndims == 2:
        // Velocity and magnetic fields
        v[0] = rcprho*u[${i}][1]; v[1] = rcprho*u[${i}][2];
        B[0] = u[${i}][3]; B[1] = u[${i}][4];

        // Compute B·v
        Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=2)};
        divB = tdivtconf[${i}][5];

        // Untransform the divergences and apply the source terms
        // Density
        tdivtconf[${i}][0] = -rcpdjac[${i}]*tdivtconf[${i}][0] + ${srcex[0]};
        // Momentum
        tdivtconf[${i}][1] = -rcpdjac[${i}]*(divB*B[0] + tdivtconf[${i}][1]) + ${srcex[1]};
        tdivtconf[${i}][2] = -rcpdjac[${i}]*(divB*B[1] + tdivtconf[${i}][2]) + ${srcex[2]};
        // Magnetic field
        tdivtconf[${i}][3] = -rcpdjac[${i}]*(divB*v[0] + tdivtconf[${i}][3]) + ${srcex[3]};
        tdivtconf[${i}][4] = -rcpdjac[${i}]*(divB*v[1] + tdivtconf[${i}][4]) + ${srcex[4]};
        // DivB
        tdivtconf[${i}][5] = 0.0; 
        // Energy
        tdivtconf[${i}][6] = -rcpdjac[${i}]*(divB*Bdotv + tdivtconf[${i}][6]) + ${srcex[6]}; 

    % elif ndims == 3:
        // Velocity and magnetic fields
        v[0] = rcprho*u[${i}][1]; v[1] = rcprho*u[${i}][2]; v[2] = rcprho*u[${i}][3];
        B[0] = u[${i}][4]; B[1] = u[${i}][5]; B[2] = u[${i}][6];

        // Compute B·v
        Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=3)};
        divB = tdivtconf[${i}][7];

        // Untransform the divergences and apply the source terms
        // Density
        tdivtconf[${i}][0] = -rcpdjac[${i}]*tdivtconf[${i}][0] + ${srcex[0]};
        // Momentum
        tdivtconf[${i}][1] = -rcpdjac[${i}]*(divB*B[0] + tdivtconf[${i}][1]) + ${srcex[1]};
        tdivtconf[${i}][2] = -rcpdjac[${i}]*(divB*B[1] + tdivtconf[${i}][2]) + ${srcex[2]};
        tdivtconf[${i}][3] = -rcpdjac[${i}]*(divB*B[2] + tdivtconf[${i}][3]) + ${srcex[3]};
        // Magnetic field
        tdivtconf[${i}][4] = -rcpdjac[${i}]*(divB*v[0] + tdivtconf[${i}][4]) + ${srcex[4]};
        tdivtconf[${i}][5] = -rcpdjac[${i}]*(divB*v[1] + tdivtconf[${i}][5]) + ${srcex[5]};
        tdivtconf[${i}][6] = -rcpdjac[${i}]*(divB*v[2] + tdivtconf[${i}][6]) + ${srcex[6]};
        // DivB
        tdivtconf[${i}][7] = 0.0; 
        // Energy
        tdivtconf[${i}][8] = -rcpdjac[${i}]*(divB*Bdotv + tdivtconf[${i}][8]) + ${srcex[8]}; 
    % endif
% endfor

// Compute forward Euler prediction of next time step
% for i,v in pyfr.ndrange(nupts, nvars):
    newsol[${i}][${v}] = u[${i}][${v}] + ${dt}*tdivtconf[${i}][${v}];
% endfor

// Get modal form at next time step
% for i,v in pyfr.ndrange(nupts, nvars):
    newsolmodes[${i}][${v}] = ${' + '.join('{jx}*newsol[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor

// Compute zeta for positivity-preserving
// ***********************************************************
    fpdtype_t zeta_low = 0;
    fpdtype_t zeta_high = 0.5;
    fpdtype_t zeta = 0;
    fpdtype_t pmin, dmin, emin, p, d, e, withinbounds;
    fpdtype_t neginf = -${large_number}; 

    % if alpha != 1.0:
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta_low', 'neginf', 'withinbounds')};
        if (withinbounds == 1){
            zeta = 0;
        }
        else {

            ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta_high', 'neginf', 'withinbounds')};
            if (withinbounds == 0){ 
                zeta_high = 1; 
                ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta_high', 'neginf', 'withinbounds')};
                if (withinbounds == 0){ 
                    zeta_high = 4; 
                    ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta_high', 'neginf', 'withinbounds')};
                    if (withinbounds == 0){ 
                        zeta_high = 10; 
                        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta_high', 'neginf', 'withinbounds')};
                        if (withinbounds == 0){ 
                            zeta_high = 50; 
                        }
                    }
                }
            }
            
            
            for (int iter = 0; iter < ${niters}; iter++) {
                zeta = 0.5*(zeta_low + zeta_high);
                ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta', 'neginf', 'withinbounds')};

                if (withinbounds == 1) {
                    zeta_high = zeta; 
                }
                else {
                    zeta_low = zeta;
                }
            }

            zeta = zeta_high;
            
        }


        // Verify that it is positivity-preserving and revert to mean mode if not
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta', 'neginf', 'withinbounds')};
        if (withinbounds == 0){
            % for i,v in pyfr.ndrange(nupts, nvars):
                filtsol_pp[${i}][${v}] = newsolmodes[0][${v}]/${mean_mode_value}; // Factor of 2 in VDM
            % endfor
        }
    % endif
// ***********************************************************

// Compute zeta for positivity-preserving and minimum entropy principle satisfying
// ***********************************************************
    zeta_low = 0;
    zeta_high = 0.5;
    zeta = 0;

    % if alpha != 0.0:
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta_low', 'ent_min', 'withinbounds')};
        if (withinbounds == 1){
            zeta = 0;
        }
        else {

            ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta_high', 'ent_min', 'withinbounds')};
            if (withinbounds == 0){ 
                zeta_high = 1; 
                ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta_high', 'ent_min', 'withinbounds')};
                if (withinbounds == 0){ 
                    zeta_high = 4; 
                    ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta_high', 'ent_min', 'withinbounds')};
                    if (withinbounds == 0){ 
                        zeta_high = 10; 
                        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta_high', 'ent_min', 'withinbounds')};
                        if (withinbounds == 0){ 
                            zeta_high = 50; 
                        }
                    }
                }
            }
            
            
            for (int iter = 0; iter < ${niters}; iter++) {
                zeta = 0.5*(zeta_low + zeta_high);
                ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta', 'ent_min', 'withinbounds')};

                if (withinbounds == 1) {
                    zeta_high = zeta; 
                }
                else {
                    zeta_low = zeta;
                }
            }

            zeta = zeta_high;
            
        }


        // Verify that it is positivity-preserving and revert to mean mode if not
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta', 'neginf', 'withinbounds')};
        if (withinbounds == 0){
            % for i,v in pyfr.ndrange(nupts, nvars):
                filtsol_mep[${i}][${v}] = newsolmodes[0][${v}]/${mean_mode_value}; // Factor of 2 in VDM
            % endfor
        }
    % endif
// ***********************************************************

// Compute convex combination (relaxed entropy principle) and forward Euler approximation of -divF
fpdtype_t sol_rep;
% for i,v in pyfr.ndrange(nupts, nvars):
    % if alpha == 0.0:
        sol_rep = filtsol_pp[${i}][${v}];
    % elif alpha == 1.0:
        sol_rep = filtsol_mep[${i}][${v}];
    % else:
        sol_rep = ${alpha}*filtsol_mep[${i}][${v}] + (1 - ${alpha})*filtsol_pp[${i}][${v}];
    % endif
    tdivtconf[${i}][${v}] = (sol_rep - u[${i}][${v}])/${dt};
% endfor



</%pyfr:kernel>

