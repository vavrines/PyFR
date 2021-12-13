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
fpdtype_t tmp;

// Compute -divF
% for i in range(nupts):
    % for v, ex in enumerate(srcex):
        tdivtconf[${i}][${v}] = -rcpdjac[${i}]*tdivtconf[${i}][${v}] + ${ex};
    % endfor
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
        
        
        % for fiter in range(niters):
            zeta = 0.5*(zeta_low + zeta_high);
            ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta', 'neginf', 'withinbounds')};

            if (withinbounds == 1) {
                zeta_high = zeta; 
            }
            else {
                zeta_low = zeta;
            }
        % endfor

        zeta = zeta_high;
        
    }


    // Verify that it is positivity-preserving and revert to mean mode if not
    ${pyfr.expand('filter', 'newsolmodes', 'filtsol_pp', 'zeta', 'neginf', 'withinbounds')};
    if (withinbounds == 0){
        % for i,v in pyfr.ndrange(nupts, nvars):
            filtsol_pp[${i}][${v}] = newsolmodes[0][${v}];
        % endfor
    }
// ***********************************************************

// Compute zeta for positivity-preserving and minimum entropy principle satisfying
// ***********************************************************
    zeta_low = 0;
    zeta_high = 0.5;
    zeta = 0;

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
        
        
        % for fiter in range(niters):
            zeta = 0.5*(zeta_low + zeta_high);
            ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta', 'ent_min', 'withinbounds')};

            if (withinbounds == 1) {
                zeta_high = zeta; 
            }
            else {
                zeta_low = zeta;
            }
        % endfor

        zeta = zeta_high;
        
    }


    // Verify that it is positivity-preserving and revert to mean mode if not
    ${pyfr.expand('filter', 'newsolmodes', 'filtsol_mep', 'zeta', 'neginf', 'withinbounds')};
    if (withinbounds == 0){
        % for i,v in pyfr.ndrange(nupts, nvars):
            filtsol_mep[${i}][${v}] = newsolmodes[0][${v}];
        % endfor
    }
// ***********************************************************

// Compute convex combination (relaxed entropy principle) and forward Euler approximation of -divF
fpdtype_t sol_rep;
% for i,v in pyfr.ndrange(nupts, nvars):
    sol_rep = ${alpha}*filtsol_mep[${i}][${v}] + (1 - ${alpha})*filtsol_pp[${i}][${v}];
    tdivtconf[${i}][${v}] = (sol_rep - u[${i}][${v}])/${dt};
% endfor

</%pyfr:kernel>


