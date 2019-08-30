function [U0,U1,U2,U3] = get_energy_dens(EV,dEV,fV,VuV,Elims)
% GET_ENERGY_DENS Get energy density of whole and parts of distribution
%
%   
%   PSD is in SI, velocity in km/s, energy in eV
%   



u = irf_units;
N = size(EV,1);

Ufun = @(E,dE,f)4*pi*sqrt(2/u.mp^3)*sum(dE.*sqrt(E.^3).*f);


% solar wind energy in [J]
EswV = .5*u.mp*sum(VuV.^2*1e6,2);

% loop over events
U0 = zeros(N,1);
U1 = zeros(N,1);
U2 = zeros(N,1);
U3 = zeros(N,1);


for ii = 1:N
    Esw = EswV(ii);
    
    % energy table in [J]
    Etemp = EV(ii,:)*u.e;
    dEtemp = dEV(ii,:)*u.e;
    % psd of event in [s^3/m^6]
    ftemp = fV(ii,:);
    
    % get total energy density
    U0(ii) = Ufun(Etemp,dEtemp,ftemp);

    % loop over limits to get partial energy densities
    for iL = 1:length(Elims)
        % find first index with upper limit over energy limit
        idElast = find(Etemp+dEtemp/2>Elims(iL)*Esw,1,'first');
        
        % check if there are at least 2 energy bins over limit
        % return 0 if not
        if idElast >= length(Etemp)-1
            continue;
        end
        
        % center energies of all bins partially or fully above limit
        Epart = Etemp(idElast:end);
        % energy bin width of all bins partially or fully above limit
        dEpart = dEtemp(idElast:end);
        
        % update first dEpart
        dEpart(1) = dEpart(1)-(Elims(iL)*Esw-Epart(1)-dEpart(1)/2);
        % get new energy center (note: linear mean)
        Epart(1) = Epart(2)-dEpart(2)/2-dEpart(1)/2;
        
        % the phase-space density is unaltered
        fpart = ftemp(idElast:end);
        
        switch iL
            case 1
                U1(ii) = Ufun(Epart,dEpart,fpart);
            case 2
                U2(ii) = Ufun(Epart,dEpart,fpart);
            case 3
                U3(ii) = Ufun(Epart,dEpart,fpart);
        end
    end
end



end

