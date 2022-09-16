# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              grad='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
    fpdtype_t inv_rho = 1/(${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))});
% for i, ex in enumerate(srcex):
% if i <= nvars - nspec:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
% else:
    tdivtconf[${i}] = -inv_rho*(${' + '.join('u[{k}]*grad[{l}][{j}]'.format(k=nspec+l, j=nspec+ndims+1+i) for l in range(ndims))}) + ${ex};
% endif
% endfor
</%pyfr:kernel>
