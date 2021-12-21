# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mhd.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.mhd.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              entminl='in view fpdtype_t'
              entmin_intl='inout view fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
fpdtype_t ur[${nvars}];
${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};
${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'artviscl', 'nl', 'magnl')};

fpdtype_t d = ur[0];
% if ndims == 2:
    rhou = ur[1]; rhov = ur[2];
    Bx = ur[3]; By = ur[4];
    E = ur[6];

    BdotB2 = 0.5*(Bx*Bx + By*By);

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov)
                             - BdotB2);
% elif ndims == 3:
    rhou = ur[1]; rhov = ur[2]; rhow = ur[3];
    Bx = ur[4]; By = ur[5]; Bz = ur[6];
    E = ur[8];

    BdotB2 = 0.5*(Bx*Bx + By*By + Bz*Bz);

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - (0.5/d)*(rhou*rhou + rhov*rhov + rhow*rhow)
                             - BdotB2);
% endif

fpdtype_t entminr = d*log(p/pow(d, ${c['gamma']}));
entmin_intl = fmin(entminl, entminr);
</%pyfr:kernel>
