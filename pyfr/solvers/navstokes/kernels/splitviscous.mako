# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='splitviscous' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf_inv='inout fpdtype_t[${str(nvars)}]'
              tdivtconf_vis='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

fpdtype_t divF_inv, divF_vis;

// Compute -divF and forward Euler prediction of next time step
// Separate viscous and inviscid components
% for v, ex in enumerate(srcex):
    divF_inv = -rcpdjac*tdivtconf_inv[${v}] + ${ex};
    divF_vis = -rcpdjac*tdivtconf_vis[${v}] + ${ex};

    // Store inviscid forward Euler prediction
     tdivtconf_inv[${v}] = u[${v}] + ${dt}*divF_inv; 

    // Store viscous difference of update in self._scal_upts_cpy
    tdivtconf_vis[${v}] = ${dt}*(divF_vis - divF_inv);
% endfor

</%pyfr:kernel>


