function fmean = sh_get_omni_dist(dist,V0)
% SH_GET_OMNI_DIST Get omnidirectional distribution in an arbitrary frame
%   
%   fOmni = SH_GET_OMNI_DIST(dist,V) returns the omnidirectional
%   distribution of the PDist object dist in the frame moving with velocity
%   V. V can either be a 1x3 vector for a frame moving with constant
%   velocity or a TSeries object for a changing frame. The unit of the
%   input distribution is conserved in the output. The unit of V is [km/s].
% 
%   In mathematical terms, this function calculates the spherical mean of
%   the distribution function https://en.wikipedia.org/wiki/Spherical_mean
%   centered on the velocity in the input. It accomplishes this with a
%   Monte-Carlo method where a number of points (5e3 per channel) evenly
%   distributed on a sphere centered on V and with radius corresponding to
%   each energy value v = sqrt(2*E/m). The omnidirectional
%   distribution is then the average value of the ditribution the
%   Monte-Carlo points for the given energy value.
%   
%
% ASSUMES IONS! ASSUMES NO ALTERNATING ENERGY TABLES! 
%
% Tested but somewhat experimental. Don't drink and use! 
% 
% See also PDIST.OMNI
% 
% Written by A. Johlander
%

% TODO:
%   - Allow alternating energy tables
%   - More inputs: energy grid, number of MC points, ...
%   - Allow for electrons
%   - Add info to the output PDist object
%   - Add to the PDist class


%% units and constants
u = irf_units;

% number of Monte-Carlo points 
nMC = 5e3; % per energy channel

%% handle inputs
% 

if nargin == 1
    V0 = [0,0,0];
end

% check species
switch dist.species
    case 'ions'
        M = u.mp;
    case 'electrons'
        M = u.me
end

nt = length(dist);

%% Handle input velocity
% make sure V0 is an array with size PDist.length
if isa(V0,'TSeries')
    %if size(V0.data,1) == size(PDist.data,1)
    if V0.length == dist.length
        irf.log('w','assuming same time sampling')
        V0v = V0.data;
    else
        irf.log('w','resampling velocity')
        V0v = V0.resample(dist).data;
    end
elseif isnumeric(V0)
    irf.log('w','using constant velocity')
    if size(V0,1)==3
        V0v = repmat(V0,nt,1);
    else
        V0v = repmat(V0,nt,1);
    end
end

% for good measure
V0v = double(V0v);

%% pre treat distribution object
% start with the normal command so it has the same structure
fmean = dist.omni;

% pre allocate data matrix
nE = size(fmean.data,2);
fmeanData = zeros(nt,nE);

% assume they are all the same
E = dist.depend{1}(1,:); % [eV]
% delta energies [eV]
dEp = dist.ancillary.delta_energy_plus(1,:);
dEm = dist.ancillary.delta_energy_minus(1,:);
% energy edges
Ee = [E-dEm,E(end)+dEp(end)]; % [eV]
ve = sqrt(2*Ee*u.e/M)*1e-3; % [km/s]


% elevation angle
th = double(dist.depend{3}); % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radians
dth = median(diff(th));
the = [th-dth/2,th(end)+dth/2]; % [radians]

% Now, energy table is the same as input distribution
Ef = E; % [eV]
nEf = nE; 
% velocity of grid in [km/s]
vf = sqrt(2*Ef*u.e/M)*1e-3; % [km/s]


%% create mc particles
% acceptence-rejection method (not so slow apparently)
thp = zeros(1,nMC);
count = 1;
while count <= nMC
    rn1 = rand(1); % [0,1] flat
    rn1 = (rn1-0.5)*pi; % [-pi/2,pi/2] flat
    
    rn2 = rand(1);
    if rn2<abs(cos(rn1))
        % accept
        thp(count) = rn1;
        count = count+1;
    end % else reject and try again
end

phip = (rand(1,nMC)-0.5)*2*pi; % [-pi,pi] flat

% velocities in spacecraft frame
vxp = zeros(nEf,nMC);
vyp = zeros(nEf,nMC);
vzp = zeros(nEf,nMC);

% use the same angle values for all energy values
for ii = 1:nEf
    [vxp(ii,:),vyp(ii,:),vzp(ii,:)] = sph2cart(phip,thp,vf(ii));
end


%% get spherical mean of distribution
% loop over time steps
fprintf('it = %4.0f/%4.0f\n',0,nt) % display progress
for ii = 1:nt
    if mod(ii,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],ii,nt); end % display progress
    
    % ---------------- instrument part ----------------
    % psd
    F3d = double(squeeze(double(dist.data(ii,:,:,:)))); % whatever units
    
    % azimuthal angle
    phi = double(dist.depend{2}(ii,:)); % in degrees
    phi = phi-180;
    phi = phi*pi/180; % in radians
    dphi = median(diff(phi));
    % edges of bins
    phie = [phi-dphi/2,phi(end)+dphi/2];
    
    % ---------------- grid part ----------------
    % velocities of mc points in desired frame
    vxp0 = vxp+V0v(ii,1); 
    vyp0 = vyp+V0v(ii,2);
    vzp0 = vzp+V0v(ii,3);
    
    % back to the (future) spherical frame
    [phip0,thp0,vp0] = cart2sph(vxp0,vyp0,vzp0);
    
    % get good indices
    idPhip = discretize(phip0,phie);
    idThp = discretize(thp0,the);
    idVp = discretize(vp0,ve);
    
    % ---------------- calculating mean psd ----------------
    for jj = 1:nEf
        % get the instrument bin indices of all MC points
        idMC = sub2ind(size(F3d),idVp(jj,:),idPhip(jj,:),idThp(jj,:));
        % sometimes there is no instrument bin corresponding to the MC
        % point, those indices become NaNs
        idMC = idMC(~isnan(idMC));
        fmeanData(ii,jj) = mean(F3d(idMC));
    end 
end

% put value in output
%fmeanData(isnan(fmeanData)) = 0;
fmean.data = fmeanData;
