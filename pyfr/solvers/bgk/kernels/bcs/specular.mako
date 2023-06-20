# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ntol = 1e-6 %>
<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>

// If left/right normal:
if (abs(abs(nl[0]) - 1.0) < ${ntol} && abs(nl[1]) < ${ntol}) {
	% for i in range(nvars):
		fr[${i}] = fl[${LRidxs[i]}];
	% endfor
}
// If up/down normal:
else if (abs(abs(nl[1]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol}) {
	% for i in range(nvars):
		fr[${i}] = fl[${UDidxs[i]}];
	% endfor
}
// If diag normal:
else if (abs(-nl[0] - nl[1]) < ${ntol}) {
	% for i in range(nvars):
		fr[${i}] = fl[${DRidxs[i]}];
	% endfor
}

</%pyfr:macro>
