function specSlope = get_spectral_slope(EV,fV,VuV,Elim)



u = irf_units;
N = size(EV,1);

% solar wind energy in [J]
EswV = .5*u.mp*sum(VuV.^2*1e6,2);

% loop over events
specSlope = zeros(N,1);

for ii = 1:N
    Esw = EswV(ii);
    % energy table for the event in [J]
    Etemp = EV(ii,:)*u.e;
    
    idE = Etemp/Esw>Elim(1) & Etemp/Esw<Elim(2);
    
    ftemp = fV(ii,idE);
    Etemp = Etemp(idE);
    
    p = polyfit(log10(Etemp),log10(ftemp),1);
    
    specSlope(ii) = p(1);
end

