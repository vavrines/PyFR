<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
% for i, ex in enumerate(srcex):
    % if i == 1 and c['omg'] != 0:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex} - ${c['omg']}*u[2];
    % elif i == 2 and c['omg'] != 0:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex} + ${c['omg']}*u[1];
    % else:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
    % endif
% endfor
</%pyfr:kernel>
