# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf_inv' ndim='2'
              t='scalar fpdtype_t'
              uout='inout fpdtype_t[${str(nvars)}]'
              uin='in fpdtype_t[${str(nvars)}]'>

% for i in range(nvars):
    uout[${i}] = (uout[${i}] - uin[${i}])/${dt};
% endfor



</%pyfr:kernel>
