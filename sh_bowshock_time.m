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
%                       (default: 'farris') (ONLY USE FARRIS NOW)
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
%       'bspos'     -   only for internal use, contains bs position  
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


% structure of bspos
%   Xsh
%   Ysh
%   Zsh
%   nSh
%   L
%   eps 
%   alpha
%   sig

%% handle input
args = varargin;
nargs = numel(varargin);

% Set default values
shModel = 'farris';
Nbs = 2e2;
Nt = 1e2;
Bu = [.5,3,-1];
Vu = [-425,30,0];
rsc = [15,0,0];
nuInp = 0;
thMax = 100;
angMin = 10; % degrees
outputSurf = 0;
debugLevel = 1;
bsposInput = 0;
bsposOutput = 0;

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
    args = [cell(1,2*Nfield),args];
    
    for ii = 1:Nfield
        args{(ii-1)*2+1} = fnms{ii};
        args{(ii-1)*2+2} = inpstr.(fnms{ii});
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
        case 'thmax'
            thMax = args{2};
        case 'angmin'
            angMin = args{2};
        case 'outputsurf'
            outputSurf = args{2};
        case 'nbs2'
            Nbs2 = args{2};
        case 'thmax2'
            thMax2 = args{2};
        case 'debug'
            debugLevel = args{2};
        case 'bspos'
            bspos = args{2};
            bsposInput = 1;
        case 'bsposout'
            bsposOutput = 1;
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


% debug header
if debugLevel >= 1
    fprintf('######################################################## \n');
    fprintf('Running sh_bowshock time... \n');
    fprintf('\n');
    fprintf('------------------------------------------------------ \n');
    fprintf(['SC position: \t [',num2str(round(rsc',2)),'] RE \n']);
    fprintf(['Upstream B: \t [',num2str(round(Bu,2)),'] nT \n']);
    fprintf(['Upstream V: \t [',num2str(round(Vu,2)),'] km/s \n']);
    fprintf('------------------------------------------------------ \n');
    fprintf(['BS model name: \t ',shModel,' \n']);
end


%% get the 3D shock model in GSE
if ~bsposInput
    [Xsh,Ysh,Zsh,nSh,L,eps,alpha,sig] = get_3D_model(rsc,shModel,Nbs,Vu,thMax);
else
    % read structure
    Xsh = bspos.Xsh;
    Ysh = bspos.Ysh;
    Zsh = bspos.Zsh;
    nSh = bspos.nSh;
    L = bspos.L;
    eps = bspos.eps;
    alpha = bspos.alpha;
    sig = bspos.sig;
end

% Write out 
if debugLevel >= 1
    
    fprintf(['\t L = \t ',num2str(L),' RE \n']);
    fprintf(['\t eps = \t ',num2str(eps),' \n']);
    fprintf(['\t alp = \t ',num2str(round(alpha,3)),' deg \n']);
    fprintf(['\t sig = \t ',num2str(sig),' \n']);
end

%% quick way to check if angle is too low compared to model
thBr = atand(sum(Bu(2:3).^2)/abs(Bu(1))); % Angle between Bu and GSE-X

if eps>1 % asymptote for hyperboloid (in natural system but whatever)
    asympAng = atand(sqrt(eps^2-1)); 
else
    asympAng = 0; %
end

% if angle is too small the starting point cannot reliably be found
if asympAng+angMin>thBr
    % do something clever if that is the case
    
    % debug
    if debugLevel >= 1
        fprintf('------------------------------------------------------ \n');
        fprintf('WARNING!!! \n'); 
        fprintf('Magnetic field is too radial to establish starting point \n'); 
        fprintf('Returning NaNs... \n'); 
        fprintf('Exiting function... \n'); 
        fprintf('\n'); 
        fprintf('\n'); 
    end
    %disp('error')
    % run again to get output structure    
    outstrct = sh_bowshock_time('debug',0,'nt',Nt);
    % special case when Tc is returned
    if outputSurf
        outstrct.X = zeros(Nbs);
        outstrct.Y = zeros(Nbs);
        outstrct.Z = zeros(Nbs);
        outstrct.Tc = zeros(Nbs);
        outstrct.cosThBn = zeros(Nbs);
    end
        
        
    fns = fieldnames(outstrct);
    for ifn = 1:length(fns)
        % put all to nans
        outstrct.(fns{ifn}) = nan(size(outstrct.(fns{ifn})));
    end
    return;
end


%% Find the spacecraft position on the bow shock
% distance between SC and BS for all points [RE]
% find closest index to sc
distBsSc = sqrt((rsc(1)-Xsh).^2+(rsc(2)-Ysh).^2+(rsc(3)-Zsh).^2);

% find the shortest distance
idmin = find(distBsSc==min(min(distBsSc)));
% idmin can be more than one for some reason, hack to fix
idmin = idmin(1);
%[id1min,id2min] = ind2sub([N,N],idmin); % to subindex
% set "new" spacecraft position, should be close to rsc
rscBS = [Xsh(idmin),Ysh(idmin),Zsh(idmin)];


%% Find starting points of fieldlines
% field lines always start at 90 degrees

% get the cosine of thBn (depends on sign of Bn)
cosThBn = zeros(Nbs,Nbs);
for ii = 1:3 % sum over components
    cosThBn = cosThBn+bu(ii)*squeeze(nSh(:,:,ii));
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


% if there are no points with 90 deg, give an error
if numel(find(id90bool))==0 
    disp('fel1')
end


%% find magnetic field/solar wind/sc plane
% vector that is perpendicular to both the magnetic field and flow velocity
out_of_plane_vector = cross(vu,bu)/norm(cross(vu,bu));

% we call the plane made up by bu and vu the "foot plane"
% get distance of bow shock to sc position in bow shock perp to the foot
% plane
distToFootPlane = out_of_plane_vector(1)*Xsh+out_of_plane_vector(2)*Ysh+...
    out_of_plane_vector(3)*Zsh-dot(out_of_plane_vector,rscBS);


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


% if there is no intersection point between 90 deg and field travel, give
% an error
if numel(find(id90bool&idFootBool))==0 
    disp('fel2')
end




%% something

% now get a vector that is perpendicular to the out-of-plane vector as well
% as the magnetic field, the fieldline will travel along this vector
field_travel_direction = cross(out_of_plane_vector,bu)/norm(cross(out_of_plane_vector,bu));

% we call the plane made up by bu and cross(bu,vu) the "field plane"
% distToFieldPlane = zeros(Nbs,Nbs);
% 
% for ii = 2:Nbs-1 % th
%     for jj = 2:Nbs-1 % phi
%         distToFieldPlane(ii,jj) = dot(field_travel_direction,[Xsh(ii,jj);Ysh(ii,jj);Zsh(ii,jj)]);
%     end
% end

distToFieldPlane = field_travel_direction(1)*Xsh+field_travel_direction(2)*Ysh+...
    field_travel_direction(3)*Zsh;


%% get total connection time (time evolution later)
% THIS IS PROBABLY WRONG!

startPointX = mean(Xsh(id90bool&idFootBool));
startPointY = mean(Ysh(id90bool&idFootBool));
startPointZ = mean(Zsh(id90bool&idFootBool));
rStart = [startPointX,startPointY,startPointZ];

% tend = dot(rsc-rStart',field_travel_direction)*u.RE/dot(Vu*1e3,field_travel_direction);

% or find distance to tangent line along Vu
disp('do not want!')


%% less lazy way
% vector pointing upstream along bu 
dhat = sign(dot(bu,vu))*bu;
% vector pointing upstream along vu 
bhat = -vu;
% field line function connecting to starting position
fieldLinePosFun = @(d)rStart+dhat*d;
% flow line function connecting to the spacecraft position
flowLinePosFun = @(b)rsc'+bhat*b;

field2flowLineDistFun = @(dbArr)norm(flowLinePosFun(dbArr(2))-fieldLinePosFun(dbArr(1)));

[dbval,~] = fminsearch(field2flowLineDistFun,[1,1]);

tend = dbval(2)*u.RE/norm(Vu*1e3);


% %% lazy way (scales with n^2)
% nn = 1e3;
% % distance along bu
% darr = linspace(0,10,nn);
% % distance along vu
% barr = linspace(0,5,nn);
% % field line points
% fieldlinepos = repmat(rStart,nn,1)'+repmat(darr,3,1).*repmat(bu,nn,1)';
% % flow line points
% flowlinepos = repmat(rsc',nn,1)'-repmat(barr,3,1).*repmat(vu,nn,1)';
% 
% minDistArr = zeros(1,nn);
% field2flowDistArr = zeros(nn,nn);
% for ii = 1:nn % b
%     field2flowDistArr(:,ii) = sqrt(sum((flowlinepos(:,ii)-fieldlinepos).^2));
%     minDistArr(ii) = min(field2flowDistArr(:,ii));
% end
% %     for ii = 1:nn % b
% %     
% %         
% %     end
% 
% % hack
if tend<0; tend = 0; end


%% get rid of "other side of field line"

if cosThBn(idmin)>=0 % thetaBn<90
    idFinal = find(idFootBool==1 & cosThBn>0);
else % thetaBn>90
    idFinal = find(idFootBool==1 & cosThBn<0);
end

% number of foot points 
%Nf = numel(find(idFootBool));

%[id1Final,id2Final] = ind2sub([N,N],idFinal);

distArray = distToFieldPlane(idFinal);
cosThBnArray = cosThBn(idFinal);
%nShArray = zeros(numel(find(idFinal)),3)
nx = nSh(:,:,1); ny = nSh(:,:,2); nz = nSh(:,:,3);
nxArray = nx(idFinal); nyArray = ny(idFinal); nzArray = nz(idFinal);
xArray = Xsh(idFinal); yArray = Ysh(idFinal); zArray = Zsh(idFinal);


%% bin values

% make nice distance vector 
% everything is flipped for now
distVec = linspace(dot(rscBS,field_travel_direction),max(distArray),Nt);
distVecDiff = median(diff(distVec));
distVecEdges = [distVec-distVecDiff/2,distVec(end)+distVecDiff/2];

rShVec = zeros(3,Nt);
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
    % give warning
    if debugLevel >= 1
        fprintf('------------------------------------------------------ \n');
        fprintf('Notice: \n'); 
        fprintf('Some parts of the vector contains NaNs. Consider increasing Nbs \n')
    end
    %disp('nans')
    cosThBnVec(idNan) = interp1(distVec(~idNan),cosThBnVec(~idNan),distVec(idNan));
    nShVec(:,idNan) = interp1(distVec(~idNan),nShVec(:,~idNan)',distVec(idNan))';
    rShVec(:,idNan) = interp1(distVec(~idNan),rShVec(:,~idNan)',distVec(idNan))';
end

% spline
cosThBnVec = smooth(cosThBnVec);

%% debug figure


if debugLevel>=2
    [h,~] = sh_figure(3,[14,16]);
    
    hca = irf_panel(h,'3dplot');
    
    % matrix for color 
    Cm = (cosThBn); %Cm(Cm>90) = 180-Cm(Cm>90);
    %Cm = Cm.*~idFootBool;
    surf(hca,Xsh,Ysh,Zsh,Cm)
    axis(hca,'equal')
    hold(hca,'on')
    shading(hca,'flat')
    hca.CLim = [-1,1];
    irf_colormap(hca,'bluered')
    
    hca.XLim = [-1,1]*max(abs(rscBS))*2+10; hca.YLim = [-1,1]*max(abs(rscBS))*2;
    hca.ZLim = [-1,1]*max(abs(rscBS))*2;
    hcb = colorbar(hca);
    hcb.Label.Interpreter = 'latex';
    hcb.Label.String = '$\cos{\theta_{Bn}}$';
    hca.View = [40,6];
    
    xlabel(hca,'X [R_E]'); ylabel(hca,'Y [R_E]'); zlabel(hca,'Z [R_E]')
    plot3(hca,rShVec(1,:),rShVec(2,:),rShVec(3,:),'r','linewidth',4)
    plot3(hca,rsc(1),rsc(2),rsc(3),'dk','markersize',8,'linewidth',3)
    plot3(hca,rscBS(1),rscBS(2),rscBS(3),'oc','markersize',8,'linewidth',3)
    
    % Bu
    arstart = rscBS+[7,0,0];
    arlen = 12;
    quiver3(hca,arstart(1),arstart(2),arstart(3),bu(1)*arlen,bu(2)*arlen,bu(3)*arlen,'linewidth',3,'MaxHeadSize',1)    
    
    sh_panel_span(h(2:end),[.1,.5])
    
    hca = irf_panel(h,'3dplot');
    plot(hca,fliplr(distVec),fliplr(cosThBnVec),'linewidth',2)
    grid(hca,'on')
    ylabel(hca,'$\cos{\theta_{Bn}}$','interpreter','latex');
    hca.FontSize = 16;
    
    hca = irf_panel(h,'3dplot');
    hold(hca,'on')
    grid(hca,'on')
    plot(hca,fliplr(distVec),(rShVec(1,:)),'k','linewidth',2)
    plot(hca,fliplr(distVec),(rShVec(2,:)),'b','linewidth',2)
    plot(hca,fliplr(distVec),(rShVec(3,:)),'r','linewidth',2)
    plot(hca,max(distVec),rscBS,'x','markersize',10,'linewidth',2)
    ylabel(hca,'$\mathbf{R}$ [$R_E$]','interpreter','latex');
    xlabel(hca,'Distance along field line path [R_E]')
    hca.FontSize = 16;
    
end



%% repeat all points on the bow shock
if debugLevel>=3 || outputSurf
    
    TconnMat = zeros(size(Xsh));
    TconnMat2 = zeros(size(Xsh));
    
    % do first point to get structure
    [Xsh2,Ysh2,Zsh2,nSh2,L2,eps2,alpha2,sig2] = get_3D_model(rsc,shModel,Nbs2,Vu,thMax2);
    % create structure
    bspos = [];
    bspos.Xsh = Xsh2;
    bspos.Ysh = Ysh2;
    bspos.Zsh = Zsh2;
    bspos.nSh = nSh2;
    bspos.L = L2;
    bspos.eps = eps2;
    bspos.alpha = alpha2;
    bspos.sig = sig2;
    
    fprintf('ii = %4.0f/%4.0f\n',0,Nbs) % display progress
    for ii = 1:Nbs % th
        fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],ii,Nbs);
        %disp([num2str(ii),'/',num2str(Nbs)])
        for jj = 1:Nbs % phi
            % whoo, recursive!
            try
                tempstrc = sh_bowshock_time('r',[Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)],'nt',Nt,'shmodel',shModel,'nbs',Nbs2,'bu',Bu,'vu',Vu,'thmax',thMax2,'debug',0,'bspos',bspos);
                TconnMat(ii,jj) = tempstrc.tend;
                TconnMat2(ii,jj) = tempstrc.t(end);
            catch
                TconnMat(ii,jj) = nan;
                TconnMat2(ii,jj) = nan;
            end
        end
    end
    
    % find 45 90 deg lines
    ln45 = nan(3,Nbs);
    thBnMat = acosd(abs(cosThBn));
    for jj = 1:Nbs-1
        for ii = 2:Nbs-1
            if (thBnMat(ii-1,jj)-45)*(thBnMat(ii+1,jj)-45)<0
                ln45(:,ii) = [Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)];
                continue;
            end
        end
    end  
elseif bsposOutput
    % do first point to get structure
    [Xsh2,Ysh2,Zsh2,nSh2,L2,eps2,alpha2,sig2] = get_3D_model(rsc,shModel,Nbs,Vu,thMax);
    % create structure
    bspos = [];
    bspos.Xsh = Xsh2;
    bspos.Ysh = Ysh2;
    bspos.Zsh = Zsh2;
    bspos.nSh = nSh2;
    bspos.L = L2;
    bspos.eps = eps2;
    bspos.alpha = alpha2;
    bspos.sig = sig2;
end

if debugLevel == 3
    disp('fooooo')
end

%% Alternative (much slower)  (far from done)
% if debugLevel>=3 || outputSurf
%     
%     yl = [-20,20];
%     zl = [-10,10];
%     
%     ny = Nbs;
%     nz = abs(diff(zl)/diff(yl))*ny;
%     
%     ye = linspace(yl(1),yl(1),ny+1);
%     ze = linspace(zl(1),zl(1),nz+1);
%     
%     dy = median(diff(ye));
%     dz = median(diff(ze));
%     
%     yc = ye(1:end-1)+dy/2;
%     zc = ze(1:end-1)+dz/2;
%     
%     
%     TconnMat = zeros(size(Xsh));
%     
%     % do first point to get structure
%     [Xsh2,Ysh2,Zsh2,nSh2,L2,eps2,alpha2,sig2] = get_3D_model(rsc,shModel,Nbs2,Vu,thMax2);
%     % create structure
%     bspos = [];
%     bspos.Xsh = Xsh2;
%     bspos.Ysh = Ysh2;
%     bspos.Zsh = Zsh2;
%     bspos.nSh = nSh2;
%     bspos.L = L2;
%     bspos.eps = eps2;
%     bspos.alpha = alpha2;
%     bspos.sig = sig2;
%     
%     fprintf('ii = %4.0f/%4.0f\n',0,Nbs) % display progress
%     for ii = 1:ny
%         fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],ii,Nbs);
%         for jj = 1:nz
%             xsh = 
%             
%             % whoo, recursive!
%             try
%                 tempstrc = sh_bowshock_time('r',[Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)],'nt',Nt,'shmodel',shModel,'nbs',Nbs2,'bu',Bu,'vu',Vu,'thmax',thMax2,'debug',0,'bspos',bspos);
%                 TconnMat(ii,jj) = tempstrc.t(end);
%             catch
%                 TconnMat(ii,jj) = nan;
%             end
%         end
%     end
%     
%     % find 45 90 deg lines
%     ln45 = nan(3,Nbs);
%     thBnMat = acosd(abs(cosThBn));
%     for jj = 1:Nbs-1
%         for ii = 2:Nbs-1
%             if (thBnMat(ii-1,jj)-45)*(thBnMat(ii+1,jj)-45)<0
%                 ln45(:,ii) = [Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)];
%                 continue;
%             end
%         end
%     end  
% end

%% derived values

% let irf_shock_parameters do everything no

% set structure
% scd = [];
% % changes that do not change
% scd.B = Bu;
% scd.V = Vu;
% if nuInp; scd.n = nu; end
% scd.Vsh = 0; % for sanity
% scd.ref_sys = 'nif';

MaVec = zeros(1,Nt);
VunVec = zeros(1,Nt);
for it = 1:Nt
    %scd.nvec = nShVec(:,it);
    %dspec = irf_shock_parameters(scd);
    
    %if isfield(dspec,'Ma'); MaVec(it) = dspec.Ma; end
    VunVec(it) = dot(nShVec(:,it),Vu);

    
end

if nuInp
    MaVec = -1e-3*(1e-9*norm(Bu)/sqrt(1e6*nu*u.mu0*u.mp))./VunVec;
end


% irf_shock_parameters does not give this for some reason
thBnVec = acosd(cosThBnVec);

% time vector in [s]
tVec = distVec*u.RE/dot(Vu*1e3,field_travel_direction);

%% set output

% remember to flip all vectors
outstrct = [];

outstrct.t = fliplr(tVec-tVec(end))';
outstrct.tend = tend;
outstrct.nvec = fliplr(nShVec)';
outstrct.thBn = flipud(thBnVec);
outstrct.Ma = fliplr(MaVec)';
outstrct.Vun = fliplr(VunVec)';
outstrct.rbs = fliplr(rShVec)';

if outputSurf
    outstrct.X = Xsh;
    outstrct.Y = Ysh;
    outstrct.Z = Zsh;
    outstrct.cosThBn = cosThBn;
    outstrct.Tc = TconnMat;
    outstrct.Tc2 = TconnMat2;
    outstrct.bspos = bspos;
end

if bsposOutput
    outstrct.bspos = bspos;
end

if debugLevel >= 1
    fprintf('------------------------------------------------------ \n');
    fprintf('The program finished without issues: \n');
    fprintf(['Total connection time: ',num2str(outstrct.t(end)),' s \n'])
    fprintf('\n')
    fprintf('\n')
end


end



%% when there is an error call this function to return empty output


%% model3d function

function [Xsh,Ysh,Zsh,nSh,L,eps,alpha,sig0] = get_3D_model(rsc,shModel,N,Vu,thMax)


%irf.log('full')
%% get model parameters (defined in irf_shock_normal)

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

r0 = [x0;y0;0];

%% get shock position, so far only in abberated system

% array of abberated angles
if thMax==0 % automatic thMax, otherwise automatic
    % elevation angle in yzx-system
    [~,thSc,~] = cart2sph(rsc(2),rsc(3),rsc(1));
    % do thMax slightly larger than spacecraft angle
    thMax = 1.5*(pi/2-thSc)*180/pi; % in degrees 
end
thpArr = linspace(-thMax,thMax,N);

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


%% get normal vectors
nSh = zeros(N,N,3);

for ii = 2:N-1 % th
    for jj = 2:N-1 % phi
        v1 = [Xsh(ii+1,jj)-Xsh(ii-1,jj),Ysh(ii+1,jj)-Ysh(ii-1,jj),Zsh(ii+1,jj)-Zsh(ii-1,jj)];
        v2 = [Xsh(ii,jj+1)-Xsh(ii,jj-1),Ysh(ii,jj+1)-Ysh(ii,jj-1),Zsh(ii,jj+1)-Zsh(ii,jj-1)];
        %%%% BUG HERE !!! ######
        %nvec = cross(v1,v2); nvec = nvec/norm(nvec); 
        % matlab function cross has terrible performance
        nvec = [v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)];
        nvec = nvec/norm(nvec);
        % point towards upstream
        %if dot(nvec,[Xsh(ii,jj),Ysh(ii,jj),Zsh(ii,jj)])<0; nvec = -nvec; end
        nvec = sign([nvec(1)*Xsh(ii,jj)+nvec(2)*Ysh(ii,jj)+nvec(3)*Zsh(ii,jj)])*nvec;
        nSh(ii,jj,:) = nvec;
        
    end
end

% hack
nSh(1,:,:) = nSh(2,:,:);
nSh(:,1,:) = nSh(:,2,:);
nSh(end,:,:) = nSh(end-1,:,:);
nSh(:,end,:) = nSh(:,end-1,:);

end


