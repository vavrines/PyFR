# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfpans' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              ku_src='inout fpdtype_t'
              wu_src='inout fpdtype_t'
              rhouforce='scalar fpdtype_t'>


% for i, ex in enumerate(srcex):
	// Turbulence sources seperately
	% if i == (nvars-2):
		tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex} + ku_src;
	% elif i == (nvars-1):
		tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex} + wu_src;
	% else:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
	% endif
% endfor

tdivtconf[${1}] += rhouforce;
tdivtconf[${nvars-3}] += rhouforce*u[${1}];


</%pyfr:kernel>
