# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<% gmrtg = (c['gamma']-1.0)/(2.0*c['gamma']) %>
<% gprtg = (c['gamma']+1.0)/(2.0*c['gamma']) %>
<% tgrgm = (2.0*c['gamma'])/(c['gamma']-1.0) %>
<% tgrgp = (2.0*c['gamma'])/(c['gamma']+1.0) %>
<% trgm = 2.0/(c['gamma']-1.0) %>
<% trgp = 2.0/(c['gamma']+1.0) %>
<% gmrgp = (c['gamma']-1.0)/(c['gamma']+1.0) %>
<% hgm = 0.5*(c['gamma']-1.0) %>
<% rgm = 1./(c['gamma']-1.0) %>
<% gamma = c['gamma'] %>

<%pyfr:macro name='lambda_phi' params='cap,p,ul,al,pl,ur,ar,pr,phi'>
    fpdtype_t fl = (p > pl) ? (p - pl)*sqrt(cap[0]/(p + cap[1])) :
                              (${trgm}*al)*(pow(p/pl,${gmrtg}) - 1.);
    fpdtype_t fr = (p > pr) ? (p - pr)*sqrt(cap[2]/(p + cap[3])) :
                              (${trgm}*ar)*(pow(p/pr,${gmrtg}) - 1.);
    phi = (fl + fr) + (ur - ul);
</%pyfr:macro>

<%pyfr:macro name='lambda_phip' params='cap,p,al,pl,ar,pr,phi'>
    fpdtype_t fl = (p > pl) ? sqrt(cap[0]/(p + cap[1]))*(1. - (p - pl)/(2.*(cap[1] + p))) :
                              (al/(${gamma}*pl))*pow(p/pl,${-gprtg});
    fpdtype_t fr = (p > pr) ? sqrt(cap[2]/(p + cap[3]))*(1. - (p - pr)/(2.*(cap[3] + p))) :
                              (ar/(${gamma}*pr))*pow(p/pr,${-gprtg});
    phi = fl + fr;
</%pyfr:macro>

<%pyfr:macro name='lambdaz' params='uz,pz,az,ps,z,v'>
     v = uz + z*az*sqrt(1. + ${gprtg}*max((ps-pz)/pz,0.));
</%pyfr:macro>

<%pyfr:macro name='lambda_newton' params='cap,ul,al,pl,ur,ar,pr,p1,p2'>
    fpdtype_t phi1,phi2,phi11,phi12,phi22,phi112,phi221;
    ${pyfr.expand('lambda_phi','cap','p1','ul','al','pl','ur','ar','pr','phi1')};
    ${pyfr.expand('lambda_phi','cap','p2','ul','al','pl','ur','ar','pr','phi2')};
    ${pyfr.expand('lambda_phip','cap','p1','al','pl','ar','pr','phi11')};
    ${pyfr.expand('lambda_phip','cap','p2','al','pl','ar','pr','phi22')};
    fpdtype_t rdp = 1./(p2 - p1);
    phi12  = (phi2  - phi1 )*rdp;
    phi112 = (phi12 - phi11)*rdp;
    phi221 = (phi22 - phi12)*rdp;
    p1 = p1 - 2.*phi1/(phi11 + sqrt(phi11*phi11 - 4.*phi1*phi112));
    p2 = p2 - 2.*phi2/(phi22 + sqrt(phi22*phi22 - 4.*phi2*phi221));
</%pyfr:macro>

<% kmax = 2 %>
<%pyfr:macro name='guermond_speed' params='rl,ul,pl,rr,ur,pr,s'>
    fpdtype_t pmin,pmax,rhomin,rhomax,cap[4],v11,v32;
    fpdtype_t rrl = 1./rl;
    fpdtype_t rrr = 1./rr;
    fpdtype_t al = sqrt(${gamma}*pl*rrl);
    fpdtype_t ar = sqrt(${gamma}*pr*rrr);
    cap[0] = ${trgp}*rrl; cap[1] = pl*${gmrgp};
    cap[2] = ${trgp}*rrr; cap[3] = pr*${gmrgp};
    if (pl < pr){
      pmin = pl; rhomin = rl;
      pmax = pr; rhomax = rr;
    }
    else{
      pmin = pr; rhomin = rr;
      pmax = pl; rhomax = rl;
    }
    fpdtype_t capAmin = ${trgp}/rhomin;
    fpdtype_t capBmin = pmin*${gmrgp};
    fpdtype_t acovmin = sqrt(${gamma}*pmin/rhomin);
    fpdtype_t acovmax = sqrt(${gamma}*pmax/rhomax);
    fpdtype_t ratio = pow(pmin/pmax,${gmrtg});
    fpdtype_t phimin = ${trgm}*acovmax*(ratio - 1.) + (ur - ul);
    if (phimin >= 0.){
        s = max(max(-(ul - al),0.), max(ur + ar,0.));
    }
    else{
        fpdtype_t phimax = (pmax-pmin)*sqrt(capAmin/(pmax + capBmin)) + (ur - ul);
        fpdtype_t ptilde = pmin*pow((acovmin + acovmax - (ur - ul)*${hgm})
                                     /(acovmin + acovmax*ratio),${tgrgm});
        fpdtype_t p1 = (phimax < 0.) ? pmax : pmin;
        fpdtype_t p2 = (phimax < 0.) ? ptilde : min(pmax,ptilde);
		% for k in range(kmax):
		        ${pyfr.expand('lambda_newton','cap','ul','al','pl','ur','ar','pr','p1','p2')};
		% endfor
        fpdtype_t pos = 1.;
        fpdtype_t neg = -1.;
        ${pyfr.expand('lambdaz','ul','pl','al','p2','neg','v11')};
        ${pyfr.expand('lambdaz','ur','pr','ar','p2','pos','v32')};
        s = max(max(v32,0.),max(-v11,0.));
    }
</%pyfr:macro>


<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rl, unl, pl, rr, unr, pr, a;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Sum the left and right velocities and take the normal
    unl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    unr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed 
    ${pyfr.expand('guermond_speed','rl','unl','pl','rr','unr','pr','a')};

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
