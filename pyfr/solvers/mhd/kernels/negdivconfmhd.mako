<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfmhd' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

% for i, ex in enumerate(srcex):    
% if i == 2*ndims + 1:
// Zero divB component of divF
tdivtconf[${i}] = 0.0;
% else:
tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
% endif
% endfor

</%pyfr:kernel>

