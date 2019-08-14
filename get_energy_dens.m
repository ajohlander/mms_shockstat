function [U0,U1,U2,U3] = get_energy_dens(EV,dEV,fV,VuV,Elims)
% GET_ENERGY_DENS Get energy density of whole and parts of distribution
%
%   
%   PSD is in SI, velocity in km/s, energy in eV
%   Can and should be improved by splitting energy bins



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
    Etemp = EV(ii,:)*u.e;
    dEtemp = dEV(ii,:)*u.e;
    ftemp = fV(ii,:);
    % get total energy density
    U0(ii) = Ufun(Etemp,dEtemp,ftemp);

    % loop over limits
    for iL = 1:length(Elims)
        % should be more sofisticated here
        Epart = Etemp(Etemp>Elims(iL)*Esw);
        dEpart = dEtemp(Etemp>Elims(iL)*Esw);
        fpart = ftemp(Etemp>Elims(iL)*Esw);
        
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

