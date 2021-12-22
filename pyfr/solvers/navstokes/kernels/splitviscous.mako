# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='splitviscous' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf_inv='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              tdivtconf_vis='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>

fpdtype_t nextsol_inv[${nupts}][${nvars}];
fpdtype_t divF_inv, divF_vis;

// Compute -divF and forward Euler prediction of next time step
// Separate viscous and inviscid components
% for i in range(nupts):
    % for v, ex in enumerate(srcex):
        divF_inv = -rcpdjac[${i}]*tdivtconf_inv[${i}][${v}] + ${ex};
        divF_vis = -rcpdjac[${i}]*tdivtconf_vis[${i}][${v}] + ${ex};

        // Store inviscid forward Euler prediction
        nextsol_inv[${i}][${v}] = u[${i}][${v}] + ${dt}*divF_inv; 

        // Store viscous component of update in self._scal_upts_cpy
        tdivtconf_vis[${i}][${v}] = ${dt}*(divF_vis - divF_inv);
    % endfor
% endfor

// Get modal form at next time step and store in self.scal_upts_outb
% for i,v in pyfr.ndrange(nupts, nvars):
    tdivtconf_inv[${i}][${v}] = ${' + '.join('{jx}*nextsol_inv[{j}][{v}]'
                                            .format(j=j, jx=jx, v=v) for j, jx in enumerate(invvdm[i]) if jx != 0)};
% endfor

</%pyfr:kernel>


