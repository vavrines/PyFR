# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>



<%pyfr:kernel name='negdivconfsplit' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf_outb='inout fpdtype_t[${str(nvars)}]'
              tdivtconf_cpy='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>


fpdtype_t tmp1, tmp2;
% for i, ex in enumerate(srcex):
    tmp1 = -rcpdjac*tdivtconf_outb[${i}] + ${ex};
    tmp2 = -rcpdjac*tdivtconf_cpy[${i}] + ${ex};

    // Swap so viscous update goes to outb
    tdivtconf_outb[${i}] = tmp2;
    tdivtconf_cpy[${i}] = tmp1;
% endfor

</%pyfr:kernel>
