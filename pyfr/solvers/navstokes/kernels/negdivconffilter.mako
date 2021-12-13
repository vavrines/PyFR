# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='filter' params='newsolmodes, filtsol, zeta, ent_min, withinbounds'>

    dmin = 999999; pmin = 999999; emin = 99999;
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

        e = d > 0 ? d*log(p/pow(d, ${c['gamma']})) : -99999;

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
fpdtype_t filtsol[${nupts}][${nvars}];
fpdtype_t filtmodes[${nupts}][${nvars}];
fpdtype_t tmp;

% for i in range(nupts):
    % for v, ex in enumerate(srcex):
        tdivtconf[${i}][${v}] = -rcpdjac[${i}]*tdivtconf[${i}][${v}] + ${ex};
    % endfor
% endfor


% for i,v in pyfr.ndrange(nupts, nvars):
    newsol[${i}][${v}] = u[${i}][${v}] + ${dt}*tdivtconf[${i}][${v}];
% endfor

% for i,v in pyfr.ndrange(nupts, nvars):
    newsolmodes[${i}][${v}] = ${' + '.join('{jx}*newsol[{j}][{v}]'.format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor


fpdtype_t zeta_low = 0;
fpdtype_t zeta_high = 0.5;
fpdtype_t zeta = 0;
fpdtype_t pmin, dmin, emin;
fpdtype_t p, d, e, withinbounds;


${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_low', 'ent_min', 'withinbounds')};
if (withinbounds == 1){
    zeta = 0;
}
else {

    ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_high', 'ent_min', 'withinbounds')};
    if (withinbounds == 0){ 
        zeta_high = 1; 
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_high', 'ent_min', 'withinbounds')};
        if (withinbounds == 0){ 
            zeta_high = 4; 
            ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_high', 'ent_min', 'withinbounds')};
            if (withinbounds == 0){ 
                zeta_high = 10; 
                ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta_high', 'ent_min', 'withinbounds')};
                if (withinbounds == 0){ 
                    zeta_high = 50; 
                }
            }
        }
    }
    
    
    % for fiter in range(niters):
        zeta = 0.5*(zeta_low + zeta_high);
        ${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta', 'ent_min', 'withinbounds')};

        if (withinbounds == 1) {
            zeta_high = zeta; 
        }
        else {
            zeta_low = zeta;
        }
    % endfor

    zeta = zeta_high;
    
}

${pyfr.expand('filter', 'newsolmodes', 'filtsol', 'zeta', 'ent_min', 'withinbounds')};

% for i,v in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${v}] = (filtsol[${i}][${v}] - u[${i}][${v}])/${dt};
% endfor
printf("%f/n", ent_min);
</%pyfr:kernel>


