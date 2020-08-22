# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='rdshocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              shockcell='out fpdtype_t'
              divf_fr='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              divf_rd='in fpdtype_t[${str(nupts)}][${str(nvars)}]'>


% if shocksensor == 'none':
	shockcell = 0;
% elif shocksensor == 'modal':
	<% se0 = crd['modal_sensor_coeff']*(order+1)**(-4*ndims) %>

	fpdtype_t totEn = 1e-15, pnEn = 1e-15, tmp;

	% for i, deg in enumerate(ubdegs):
		tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=svar)
							for j, jx in enumerate(invvdm[i]) if jx != 0)};
		totEn += tmp*tmp;
		% if deg >= order:
			pnEn += tmp*tmp;
		% endif
	% endfor

    fpdtype_t se  = pnEn/totEn;
    shockcell = (se < ${se0}) ? 0 : 1;
% else:
	shockcell = 1;
% endif

if (shockcell == 1) {
	% for i in range(nupts):
		% for j in range(nvars):
			divf_fr[${i}][${j}] = divf_rd[${i}][${j}];
		% endfor
	% endfor
}

</%pyfr:kernel>