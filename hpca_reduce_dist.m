function [dist,moms] = hpca_reduce_dist(hpca_dist,energy,azimuth,elevation,varargin)
%HPCA_REDUCE_DIST Reduce distribution function from MMS-HPCA
%
%   WORK IN PROGRESS! LOTS OF UNDOCUMENTED FEATURES AND POSSIBLE BUGS!
%
%   [cartdist,moms] = HPCA_REDUCE_DIST(DIST,E,AZ,ELE,'Opt1',OptVal1,...)
%   Rebins HPCA ion distribution DIST to a cartesian system. Also calculates
%   reduced 2D and 1D distribution functions and derived moments. E is
%   energy table, AZ is a TSeries object with azimuthal angels, and ELE is
%   a vector with elevation angles.
%       
%   Structure of output:
%      cartdist:
%           v       -   velocity in [km/s]
%           t       -   time stamps (start time)
%           xhat    -   normal vector of x direction
%           yhat    -   normal vector of y direction
%           zhat    -   normal vector of z direction
%           f       -   3D cartesian phase-space density [SI]
%           fx      -   reduced phase-space density along xhat [SI]
%           fy      -   reduced phase-space density along yhat [SI]
%           fz      -   reduced phase-space density along zhat [SI]
%           fxy     -   2D reduced phase-space density in x-y plane [SI]
%           fxz     -   2D reduced phase-space density in x-z plane [SI]
%           fyz     -   2D reduced phase-space density in y-z plane [SI]
%   
%      moms:
%           n       -   ion number density [cm^-3]
%           V       -   ion flow velocity [km/s]
%
%   Options for output:
%       'xyz'       -   3x3 matrix with [x;y;z]. z is normal to the plotted plane and
%                   x and y are made orthogonal to z and each other if they are 
%                    not already. If you want to plot different planes you have to
%                   rotate this matrix -> [y;z;x] -> [z;x;y]. Make sure the system
%                   is right-handed.
%       'nMC'       -   number of Monte Carlo iterations used for integration,
%                   default is 1000
%       'vg'        -   array with center values for the projection velocity
%                   grid in [km/s]. Some default value.
%       'm'         -   mass number of ion species.
%

%% Get velocity

u = irf_units;



%% set input

% default values
% mc number
nMC = 1000;
% new binning
dvg = 50e3; % m/s
% centers
%vpx = -max(v)*1.1:dvp:max(v)*1.1;
vgx = -2.7e6*1.1:dvg:2.7e6;
% outputmode?
% mass?
M = 1;
xphat = [1,0,0];
yphat = [0,1,0];
zphat = [0,0,1];

% check for flags
args = varargin;
nargs = length(varargin);

% loop to check for flags
have_options = nargs > 1;
while have_options
    switch(lower(args{1}))
        case 'vg'
            vgx = args{2}*1e3; % m/s
            dvg = median(diff(vgx));
        case 'nmc'
            nMC = args{2};
        %case 'dim'
         %   distDim = args{2};
        case 'xyz' % 3x3 matrix where xphat = coord_sys(1,:),...
            coord_sys = args{2};
            xphat = coord_sys(1,:)/norm(coord_sys(1,:)); % new x-axis (x-primed)
            zphat = coord_sys(3,:)/norm(coord_sys(3,:)); % vector to integrate along (z-primed)
            % only x and z are relvant
            yphat = cross(zphat,xphat); yphat = yphat/norm(yphat);
            if abs(acosd(yphat*(coord_sys(2,:)/norm(coord_sys(2,:)))'))>1
                irf.log('warning',['y (perp1) changed from [' num2str(coord_sys(2,:)/norm(coord_sys(2,:)),'% .2f') '] to [' num2str(yphat,'% .2f') '].']);
            end
        case 'm'
            M = args{2};
    end
    args = args(3:end);
    if isempty(args), break, end
end




%% 


v = sqrt(2*energy*u.e/u.mp/M); % si
% edges
ve = [v(1)-(v(2)-v(1))/2;v(2:end)-diff(v);v(end)+(v(end)-v(end-1))/2];


% edges
vpxe = [vgx-dvg/2,vgx(end)+dvg];
% y and z are the same
Np = length(vgx);

nt = length(azimuth.time);
dens = zeros(nt,1);
vel = zeros(nt,3);
F1Px = zeros(nt,Np);
F1Py = zeros(nt,Np);
F1Pz = zeros(nt,Np);
F2Pxy = zeros(nt,Np,Np);
F2Pxz = zeros(nt,Np,Np);
F2Pyz = zeros(nt,Np,Np);
F3P = zeros(nt,Np,Np,Np);

P = zeros(nt,16,16,63);


for idt = 1:nt-1
    
    it = find(hpca_dist.time.epochUnix<=azimuth.time(idt).epochUnix+0.0001 &...
        hpca_dist.time.epochUnix>=azimuth.time(idt).epochUnix-0.0001);
    
    disp([num2str(idt),'/',num2str(nt)])
    
    % phase-space density [cm^-6 s^3]
    f = hpca_dist.data(it:it+15,:,:,:);
    
    % convert to [m^-6 s^3]
    f = f*1e12;
    
    % think it's start time
    phi = squeeze(azimuth.data(idt,:,:))*pi/180;
    
    th = 90-elevation;
    th = th*pi/180;
    
    nV = length(v);
    nPhi = length(phi);
    nTh = length(th);
    
    Nbins = nV*nPhi*nTh;
    
    
    dPhi = (11.25+7)/180*pi;
    dTh = 24/180*pi;
    %dPhi = 2*pi/nPhi;
    %dTh = 2*pi/nTh;
    dV = diff(ve);
    
    vtev = zeros(nV,2);
    for iv = 1:nV
        vtev(iv,:) = [ve(iv),ve(iv+1)];
    end
    
    % % [phi,th,E] (right?)
    V = repmat(v,1,nPhi,nTh);
    V = permute(V,[2,3,1]);
    PHI = repmat(phi,1,1,nV);
    TH = repmat(th,nPhi,1,nV);
    
    % super large arrays holding all mc particles
    VX = zeros(1,Nbins*nMC);
    VY = zeros(1,Nbins*nMC);
    VZ = zeros(1,Nbins*nMC);
    
    F = zeros(1,Nbins*nMC); % psd
    
    id = 1;
    for jj = 1:nPhi % actually time steps
        for kk = 1:nTh
            for ll = 1:nV
                
                % generate random particles
                % first is not random
                dV_MC = [0;(rand(nMC-1,1)-.5)*dV(ll)];
                dPHI_MC = [0;(rand(nMC-1,1)-.5)*dPhi];
                dTH_MC = [0;(rand(nMC-1,1)-.5)*dTh];
                
                % convert instrument bin to cartesian velocity
                [vx,vy,vz] = sph2cart(PHI(jj,kk,ll)+dPHI_MC,TH(jj,kk,ll)+dTH_MC,V(jj,kk,ll)+dV_MC);
                % bit of a hack to find index
                %id = find(VX==0,1,'first');
                
                %if isempty(id); id = 1; end
                VX(id:id+nMC-1) = [vx,vy,vz]*xphat';
                VY(id:id+nMC-1) = [vx,vy,vz]*yphat';
                VZ(id:id+nMC-1) = [vx,vy,vz]*zphat';
                
                F(id:id+nMC-1) = f(jj,kk,ll);
                
                % update id 
                id = id+nMC;
            end
        end
    end
    
    
    %% Awesome cartesian rebinning

    iVx = discretize(VX,vpxe);
    iVy = discretize(VY,vpxe);
    iVz = discretize(VZ,vpxe);
    
    fp = zeros(Np,Np,Np); % average
    npp = zeros(Np,Np,Np);
    
    
    IND = sub2ind(size(fp),iVx,iVy,iVz);
    
    for ii = 1:length(IND)
        if ~isnan(IND(ii))
            fp(IND(ii)) = fp(IND(ii))+F(ii);
            npp(IND(ii)) = npp(IND(ii))+1;
        end
    end
    fp(npp~=0) = fp(npp~=0)./npp(npp~=0);
    
    
    %% reduce distribution if requested
    
    fp2dxy = sum(fp*dvg,3);
    fp2dxz = squeeze(sum(fp*dvg,2));
    fp2dyz = squeeze(sum(fp*dvg,1));
    fp1dx = sum(sum(fp*dvg^2,2),3);
    fp1dy = sum(sum(fp*dvg^2,1),3);
    fp1dz = squeeze(sum(sum(fp*dvg^2,1),2));
    
    
    
    %% calculate moments
    dens(idt) = sum(sum(sum(fp*dvg^3)));
    
    for ii = 1:Np
        for jj = 1:Np
            for kk = 1:Np
                vel(idt,:) = vel(idt,:)+[vgx(ii),vgx(jj),vgx(kk)]*fp(ii,jj,kk)*dvg^3;
            end
        end
    end
    
    vel(idt,:) = vel(idt,:)/dens(idt);
    
    F1Px(idt,:) = fp1dx;
    F1Py(idt,:) = fp1dy;
    F1Pz(idt,:) = fp1dz;
    F2Pxy(idt,:,:) = fp2dxy;
    F2Pxz(idt,:,:) = fp2dxz;
    F2Pyz(idt,:,:) = fp2dyz;
    F3P(idt,:,:,:) = fp;
    
    
end

dens(end) = nan;
n = irf.ts_scalar(azimuth.time,dens*1e-6);
v = irf.ts_vec_xyz(azimuth.time,vel*1e-3);

moms = [];
moms.n = n;
moms.V = v;

dist = [];
dist.v = vgx*1e-3;
dist.t = azimuth.time;

dist.vec1 = xphat;
dist.vec2 = yphat;
dist.vec3 = zphat;

dist.fx = F1Px;
dist.fy = F1Py;
dist.fz = F1Pz;
dist.fxy = F2Pxy;
dist.fxz = F2Pxz;
dist.fyz = F2Pyz;
dist.f = F3P;

dist.P = P;

end

