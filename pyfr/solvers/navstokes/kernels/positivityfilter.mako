# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% large_number = 10**10 %>

<%pyfr:macro name='filter_pp' params='u_modes, filtsol, zeta, withinbounds'>

    dmin = ${large_number}; pmin = ${large_number};;
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

        dmin = fmin(dmin, d);
        pmin = fmin(pmin, p);
    % endfor

    if (dmin >= ${dtol} && pmin >= ${ptol}) {
        withinbounds = 1; 
    }
    else {
        withinbounds = 0;
    }
</%pyfr:macro>


<%pyfr:kernel name='positivityfilter' ndim='1'
              t='scalar fpdtype_t'
              u_inv='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              du_vis='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>


fpdtype_t u_modes[${nupts}][${nvars}];

fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];

fpdtype_t zeta_low = 0;
fpdtype_t zeta_high = 12;
fpdtype_t zeta = 0;
fpdtype_t pmin, dmin, p, d, withinbounds;

// Add viscous component and check if positivity-preserving
% for i,v in pyfr.ndrange(nupts, nvars):
    u_inv[${i}][${v}] += du_vis[${i}][${v}];
% endfor

// Get modal form of filtered viscous solution
% for i,v in pyfr.ndrange(nupts, nvars):
    u_modes[${i}][${v}] = ${' + '.join('{jx}*u_inv[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor

${pyfr.expand('filter_pp', 'u_modes', 'u_inv', 'zeta', 'withinbounds')};

// If not positivity-preserving, continue increasing filter strength
if (withinbounds == 0) {

    for (int iter = 0; iter < ${niters}; iter++) {
        zeta = 0.5*(zeta_low + zeta_high);
        ${pyfr.expand('filter_pp', 'u_modes', 'u_inv', 'zeta', 'withinbounds')};

        if (withinbounds == 1) {
            zeta_high = zeta; 
        }
        else {
            zeta_low = zeta;
        }
    }
    zeta = zeta_high;
    
    ${pyfr.expand('filter_pp', 'u_modes', 'u_inv', 'zeta', 'withinbounds')};
}


</%pyfr:kernel>


