# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
              
fpdtype_t ur[${nvars}];
${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};
${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'artviscl', 'nl', 'magnl')};

fpdtype_t artviscr;
fpdtype_t d = ur[0];
% if ndims == 2:
    fpdtype_t p = ${c['gamma'] - 1}*(ur[${nvars - 1}] - 
                            (0.5/d)*(ur[1]*ur[1] + 
                                     ur[2]*ur[2]));
% elif ndims == 3:
    fpdtype_t p = ${c['gamma'] - 1}*(ur[${nvars - 1}] - 
                            (0.5/d)*(ur[1]*ur[1] + 
                                     ur[2]*ur[2] + 
                                     ur[3]*ur[3]));
% endif

fpdtype_t artviscr = d*log(p/pow(d, ${c['gamma']}));
artviscl = fmin(artviscl, artviscr);
</%pyfr:kernel>
