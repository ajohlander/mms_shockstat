function [outstrct] = sh_bowshock_time(varargin)
% SH_BOWSHOCK_TIME Get the history of a fieldline on the bow shock.
%   
%   WARNING: The function makes assumptions that the solar wind parameters
%   and bow shock shape are constant over time, typically a few minutes.
%   Make sure this constraint is fulfilled.
%
%   bsspec = SH_BOWSHOCK_TIME(flag1,value1,flag2,value2,...) returns a
%   structure bsspec containing the time evolution of the foot point of a
%   magnetic field line connected to the bow shock before it reaches a
%   spacecraft. 
%   
%   The fields in the output bsspec are:
%       t       -   time vector in [s], t=0 corresponds to the first
%                   encounter of the fieldline to the bow shock (thBn=90
%                   degrees at t=0)
%       nvec    -   matrix of normal unit vectors as a function of time
%       thBn    -   shock angle in [deg]
%       Ma      -   Alfven Mach number (zero if nu not given in input)
%       Mms     -   fast magnetosonic Mach number NOT IMPLEMENTED
%       Vun     -   upstream speed along normal vector
%       rbs     -   position of the field line foot point in [RE]
% 
%   The input flags are specified below, all parameters have default values
%   but it is recommended to specify all of them (not case sensitive):
%       'shModel'   -   which model to use, see irf_shock_normal for a list
%                       (default: 'farris')
%       'rsc'       -   position of the spacecraft in [RE]
%       'Bu'        -   upstream magnetic field vector in [nT]
%       'Vu'        -   upstream bulk velocity vector in [km/s]
%       'nu'        -   upstream plasma denisty in [cm^-3], needed for Mach
%                       numbers
%       'Tiu'       -   upstream ion temperature in [eV] NOT IMPLEMENTED
%       'Teu'       -   upstream electron temperature in [eV] NOT IMPLEMENTED
%       'Nbs'       -   number of points used to calculate bow shock
%                       position and normal, increase to increase accuracy,
%                       note that performance goes as Nbs^-3 (default: 200)
%       'Nt'        -   how many time steps the output is given in
%                       (default: 100)
%       'debug'     -   boolean value to display a debug figure NOT
%                       IMPLEMENTED
%
%   SH_BOWSHOCK_TIME(scd,flag1,value1,...) parts or all of the input can be
%   a structure with the flags as field names, makes function compatible
%   with irf_shock_normal and irf_shock_parameters
% 
%   
%   See also: IRF_SHOCK_NORMAL, IRF_SHOCK_PARAMETERS
%
%
%   TODO: fast Mach number, debug figure, information logs, cleanup,
%   testing



%% handle input
args = varargin;
nargs = numel(varargin);

% Set default values
shModel = 'farris';
Nbs = 2e2;
Nt = 1e2;
Bu = [3,3,0];
Vu = [-425,30,0];
rsc = [15,0,0];
nuInp = 0;


% check if input is a structure
if nargin>0 && isstruct(args{1})
    % reset args
    nargs = nargs-1;
    inpstr = args{1};
    if nargs>1; args = args(2:end); else args = {}; end
    
    
    % fieldnames
    fnms = fieldnames(inpstr);
    % number of fields
    Nfield = numel(fnms);
    % lengthen args
    args = [args,cell(1,2*Nfield)];
    
    for ii = 1:Nfield
        args{nargs+(ii-1)*2+1} = fnms{ii};
        args{nargs+(ii-1)*2+2} = inpstr.(fnms{ii});
    end
    
    nargs = nargs+Nfield;
end
    
% handle flags
have_options = nargs > 1;
while have_options
    switch(lower(args{1}))
        case 'shmodel'
            shModel = args{2};
        case 'nbs'
            Nbs = args{2};
        case 'nt'
            Nt = args{2};
        case {'b','bu'}
            Bu = args{2};
        case {'v','vu'}
            Vu = args{2};
        case {'rsc','r'}
            rsc = args{2};
        case 'nu'
            nuInp = 1;
            nu = args{2};
    end
    args = args(3:end);
    if isempty(args), break, end
end


u = irf_units;


% sc position in GSE in Earth radii must be colum vector
if size(rsc,1)<size(rsc,2)
    rsc = rsc';
end


% unit vectors of upstream B and V
bu = Bu/norm(Bu);
vu = Vu/norm(Vu); 


%% get the 3D shock model in GSE
[Xsh,Ysh,Zsh,nSh] = get_3D_model(rsc,shModel,Nbs,Vu);


%% 
% distance between SC and BS for all points [RE]
distBsSc = zeros(Nbs,Nbs);
% find closest index to sc
for ii = 1:Nbs % th
    for jj = 1:Nbs % phi
        distBsSc(ii,jj) = norm(rsc-[Xsh(ii,jj);Ysh(ii,jj);Zsh(ii,jj)]); 
    end
end

% find the shortest distance
idmin = find(distBsSc==min(min(distBsSc)));
% idmin can be more than one for some reason, hack to fix
idmin = idmin(1);
%[id1min,id2min] = ind2sub([N,N],idmin); % to subindex
% set "new" spacecraft position, should be close to rsc
rscBS = [Xsh(idmin),Ysh(idmin),Zsh(idmin)];


%% Find starting points of fieldlines
% field lines always start at 90 degrees

% get the cosine of thBn
cosThBn = zeros(Nbs,Nbs);
for ii = 1:Nbs % th
    for jj = 1:Nbs % phi
        cosThBn(ii,jj) = dot(bu,squeeze(nSh(ii,jj,:)));
    end
end


% find all points where thBn is 90 degrees
id90bool = zeros(Nbs,Nbs);
for ii = 2:Nbs-1 % th
    for jj = 2:Nbs-1 % phi
        if cosThBn(ii-1,jj)*cosThBn(ii+1,jj)<=0 || cosThBn(ii,jj-1)*cosThBn(ii,jj+1)<=0
            id90bool(ii,jj) = 1;
        end
    end
end



%% find magnetic field/solar wind/sc plane
% vector that is perpendicular to both the magnetic field and flow velocity
out_of_plane_vector = cross(vu,bu)/norm(cross(vu,bu));

% we call the plane made up by bu and vu the "foot plane"
distToFootPlane = zeros(Nbs,Nbs);

% get distance of bow shock to sc position in bow shock perp to the foot
% plane
for ii = 1:Nbs % th
    for jj = 1:Nbs % phi
        distToFootPlane(ii,jj) = dot(out_of_plane_vector,[Xsh(ii,jj);Ysh(ii,jj);Zsh(ii,jj)])-dot(out_of_plane_vector,rscBS);
    end
end


% find where this distance changes sign, these are the foot points of the
% field line 
idFootBool = zeros(Nbs,Nbs);
for ii = 2:Nbs-1 % th
    for jj = 2:Nbs-1 % phi
        if distToFootPlane(ii-1,jj)*distToFootPlane(ii+1,jj)<=0 || distToFootPlane(ii,jj-1)*distToFootPlane(ii,jj+1)<=0
            idFootBool(ii,jj) = 1;
        end
    end
end


%% something

% now get a vector that is perpendicular to the out-of-plane vector as well
% as the magnetic field, the fieldline will travel along this vector
avec = cross(out_of_plane_vector,bu)/norm(cross(out_of_plane_vector,bu));

% we call the plane made up by bu and cross(bu,vu) the "field plane"
distToFieldPlane = zeros(Nbs,Nbs);

for ii = 2:Nbs-1 % th
    for jj = 2:Nbs-1 % phi
        distToFieldPlane(ii,jj) = dot(avec,[Xsh(ii,jj);Ysh(ii,jj);Zsh(ii,jj)]);
    end
end

%% get subscript indices to foot points

% % number of foot points 
% Nf = numel(find(idFootBool));
% [id1foot,id2foot] = ind2sub([N,N],find(idFootBool));
% 
% Xf = Xsh(idFootBool==1);
% Yf = Ysh(idFootBool==1);
% Zf = Zsh(idFootBool==1);
% 
% 
% %% plot
% figure;
% foo = zeros(N,N); foo(idmin) = 1;
% surf(Xsh,Ysh,Zsh,idFootBool-id90bool+foo)
% hold on;
% plot3(rsc(1),rsc(2),rsc(3),'r*','markersize',20)
% plot3([0,20],[0,0],[0,0],'k')
% 
% axis equal;
% %colorbar
% xlim([-20,20])
% ylim([-30,30])
% zlim([-30,30])
% shading flat
% 
% %
% xlabel('X [R_E]')
% ylabel('Y [R_E]')
% zlabel('Z [R_E]')
% 
% %%
% figure
% scatter(Xsh(idFootBool==1),acosd(cosThBn(idFootBool==1)))
% hold on
% scatter(Xsh(find(foo)),acosd(cosThBn(find(foo))),'x')
% scatter(Xsh(find(idFootBool&id90bool)),acosd(cosThBn(idFootBool&id90bool)),'x')
% 
% %%
% figure
% scatter(distToFieldPlane(find(idFootBool)),acosd(cosThBn(find(idFootBool))))
% hold on
% scatter(distToFieldPlane(find(foo)),acosd(cosThBn(find(foo))),'x')
% scatter(distToFieldPlane(find(idFootBool&id90bool)),acosd(cosThBn(idFootBool&id90bool)),'x')
% plot(dot(rsc,avec)*[1,1],180*[0,1])


%% get rid of "other side of field line"

if cosThBn(idmin)>=0 % thetaBn<90
    idFinal = find(idFootBool==1 & cosThBn>0);
else % thetaBn>90
    idFinal = find(idFootBool==1 & cosThBn<0);
end

% number of foot points 
Nf = numel(find(idFootBool));

%[id1Final,id2Final] = ind2sub([N,N],idFinal);

distArray = distToFieldPlane(idFinal);
cosThBnArray = cosThBn(idFinal);
%nShArray = zeros(numel(find(idFinal)),3)
nx = nSh(:,:,1); ny = nSh(:,:,2); nz = nSh(:,:,3);
nxArray = nx(idFinal); nyArray = ny(idFinal); nzArray = nz(idFinal);
xArray = Xsh(idFinal); yArray = Ysh(idFinal); zArray = Zsh(idFinal);
disp('in a garden in the house of love');


%% bin values

% make nice distance vector 
% everything is flipped for now
distVec = linspace(dot(rscBS,avec),max(distArray),Nt);
distVecDiff = median(diff(distVec));
distVecEdges = [distVec-distVecDiff/2,distVec(end)+distVecDiff/2];

cosThBnVec = zeros(1,Nt);
nShVec = zeros(3,Nt);
% discretize the data into bins
idDistBins = discretize(distArray,distVecEdges);

for ii = 1:Nt
    nShVec(:,ii) = [mean(nxArray(idDistBins==ii)),mean(nyArray(idDistBins==ii)),mean(nzArray(idDistBins==ii))];
    rShVec(:,ii) = [mean(xArray(idDistBins==ii)),mean(yArray(idDistBins==ii)),mean(zArray(idDistBins==ii))];
    cosThBnVec(ii) = mean(cosThBnArray(idDistBins==ii));
end

% fix nans
idNan = isnan(cosThBnVec);
if ~isempty(find(idNan,1))
    irf.log('w','Some parts of the vector contains NaNs. Consider improving resolution')
    cosThBnVec(idNan) = interp1(distVec(~idNan),cosThBnVec(~idNan),distVec(idNan));
    nShVec(:,idNan) = interp1(distVec(~idNan),nShVec(:,~idNan)',distVec(idNan))';
    rShVec(:,idNan) = interp1(distVec(~idNan),rShVec(:,~idNan)',distVec(idNan))';
end

% spline
cosThBnVec = smooth(cosThBnVec);

%% derived values
disp('in love with a psycho')

% let irf_shock_parameters do everything

% set structure
scd = [];
% changes that do not change
scd.B = Bu;
scd.V = Vu;
if nuInp; scd.n = nu; end
scd.Vsh = 0; % for sanity
scd.ref_sys = 'nif';

MaVec = zeros(1,Nt);
VunVec = zeros(1,Nt);
for it = 1:Nt
    scd.nvec = nShVec(:,it);
    dspec = irf_shock_parameters(scd);
    
    if isfield(dspec,'Ma'); MaVec(it) = dspec.Ma; end
    VunVec(it) = dot(scd.nvec,Vu);
end
% irf_shock_parameters does not give this for some reason
thBnVec = acosd(cosThBnVec);

% time vector in [s]
u = irf_units;
tVec = distVec*u.RE/dot(Vu*1e3,avec);

%% set output

% remember to flip all vectors
outstrct = [];

outstrct.t = fliplr(tVec-tVec(end))';
outstrct.nvec = fliplr(nShVec)';
outstrct.thBn = flipud(thBnVec);
outstrct.Ma = fliplr(MaVec)';
outstrct.Vun = fliplr(VunVec)';
outstrct.rbs = fliplr(rShVec)';

disp('im only trying to remind you')

end

%% model3d function

function [Xsh,Ysh,Zsh,nSh] = get_3D_model(rsc,shModel,N,Vu)


%irf.log('full')
%% get model parameters (defined in irf_shock_normal)
z0 = 0; % always?

% get sigma and everything else
scd = []; scd.Bu = [0,0,0]; scd.Vu = Vu; scd.nu = 1;
scd.Bd = [0,0,0]; scd.Vd = [0,0,0]; scd.nd = 1; scd.R = rsc;

% call here to avoid double definition
dspec = irf_shock_normal(scd);

sig0 = dspec.info.sig.(shModel);
alpha = dspec.info.alpha.(shModel);
eps = dspec.info.eps.(shModel);
L = dspec.info.L.(shModel);
x0 = dspec.info.x0.(shModel);
y0 = dspec.info.y0.(shModel);

r0 = [x0;y0;z0];


%% get shock position, so far only in abberated system

% array of abberated angles
% the limits are emprical, for elliptic conic section it could be anything
thpArr = linspace(-100,100,N);

% array of radii
rpArr = sig0*L./(1+eps*cosd(thpArr));

for ii = 1:N % should be vectorized
    rpv = [rpArr(ii)*cosd(thpArr(ii));rpArr(ii)*sind(thpArr(ii))];
    %rv = (rpv+sig0*r0);
    %thv(ii) = atand(rv(2)/rv(1));
    xpsh(ii) = rpv(1); ypsh(ii) = rpv(2);
    %[xsh(ii),ysh(ii)] = pol2cart(thv(ii)
end


%% get 3d bs position in abberated system (only this!)
% rotation angles (centers?)
phiArr = linspace(0,2*pi-2*pi/N,N);
% array of abberated angles (centers)
% the limits are emprical, for elliptic conic section it could be anything


% array of radii (2D)
rpArr = sig0*L./(1+eps*cosd(thpArr));

% initilize matrices
Xpsh = zeros(N,N);
Ypsh = zeros(N,N);
Zpsh = zeros(N,N);

for ii = 1:N % th
    for jj = 1:N % phi
        % 2d
        rpv = [rpArr(ii)*cosd(thpArr(ii));rpArr(ii)*sind(thpArr(ii))];
        %rv = rpv+sig0*r0;
        %thv(ii) = atand(rv(2)/rv(1));
        Xpsh(ii,jj) = rpv(1); Ypsh(ii,jj) = rpv(2)*cos(phiArr(jj)); Zpsh(ii,jj) = rpv(2)*sin(phiArr(jj));
    end
end


% keep
%Xsh = Xpsh; Ysh = Ypsh; Zsh = Zpsh;
%R = [cosd(alpha),-sind(alpha),0;sind(alpha),cosd(alpha),0;0,0,1];


%% transform to real system
% unit vector pointing from origin of abberated system to nose of bow shock
%ehat = [sqrt(1-tand(alpha)^2/(1+tand(alpha)^2)),-sqrt(tand(alpha)^2/(1+tand(alpha)^2)),0];

% rotation matrix (as in 10.23)
R = [cosd(alpha),-sind(alpha),0;sind(alpha),cosd(alpha),0;0,0,1];
%R = [cosd(alpha),-sind(alpha);sind(alpha),cosd(alpha)];

% initilize matrices
Xsh = zeros(N,N);
Ysh = zeros(N,N);
Zsh = zeros(N,N);

for ii = 1:N % th
    for jj = 1:N % phi
        Rnat = R\[Xpsh(ii,jj);Ypsh(ii,jj);Zpsh(ii,jj)]+sig0*r0;
        
        Xsh(ii,jj) = Rnat(1); Ysh(ii,jj) = Rnat(2); Zsh(ii,jj) = Rnat(3);
        %[Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)] 
    end
end

disp('pulsar')

%% get normal vectors
nSh = zeros(N,N,3);

for ii = 2:N-1 % th
    for jj = 2:N-1 % phi
        v1 = [Xsh(ii+1,jj)-Xsh(ii-1,jj),Ysh(ii+1,jj)-Ysh(ii-1,jj),Zsh(ii+1,jj)-Zsh(ii-1,jj)];
        v2 = [Xsh(ii,jj+1)-Xsh(ii,jj-1),Ysh(ii,jj+1)-Ysh(ii,jj-1),Zsh(ii,jj+1)-Zsh(ii,jj-1)];
        
        nvec = cross(v1,v2)/norm(cross(v1,v2));
        % point towards upstream
        if dot(nvec,[Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)])<0; nvec = -nvec; end
        
        nSh(ii,jj,:) = nvec;
        
    end
end

% hack
nSh(1,:,:) = nSh(2,:,:);
nSh(:,1,:) = nSh(:,2,:);
nSh(end,:,:) = nSh(end-1,:,:);
nSh(:,end,:) = nSh(:,end-1,:);

end


