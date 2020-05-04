# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfav' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              artvisc='in fpdtype_t[${str(nvars)}]'>
% for i, ex in enumerate(srcex):
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
% endfor

% for i in range(nvars):
    tdivtconf[${i}] += artvisc[${i}];
    //if (abs(artvisc[${i}]) > 0.0000001){
    //	printf("%f\n", artvisc[${i}]);
    //}

% endfor
</%pyfr:kernel>
