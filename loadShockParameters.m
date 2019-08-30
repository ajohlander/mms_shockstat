
%% Select data file
dataFileArr = dir('*.mat');

if isempty(dataFileArr)
    irf.log('c','No data files')
    
else
    preSelectedFile = 1;
    fprintf('Data files:\n')
    for ii = 1:length(dataFileArr)
        fprintf([num2str(ii),':     ',dataFileArr(ii).name,'\n'])
        % make shockParams default file
        if strcmp(dataFileArr(ii).name,'shockParams.mat')
            preSelectedFile = ii;
        end
    end
    fprintf('\n')
    dataFileInd = irf_ask('Select data file: [%]>','listFileName',preSelectedFile);
    load(dataFileArr(dataFileInd).name)
    fprintf(['Loaded data from file: ',dataFileArr(dataFileInd).name,' \n\n'])
end



%% Select which data
%instrumentMode = irf_ask('Use EIS or only FPI? (1:FPI, 2:EIS): [%]>','instrumentMode',1);
fprintf('Put a limit on the relative energy range of the resolution: \n')
EmaxFac = irf_ask('Emax must be at least this multiple greater than Esw: [%]>','EmaxFac',20);

% select which shock model to use (all are saved)
shModel = {'farris','slho','per','fa4o','fan4o','foun'};

fprintf('\n')
fprintf('Select which bow shock model to use (only "farris" works properly)\n')
for jj = 1:length(shModel); fprintf([num2str(jj),': ',shModel{jj},'\n']); end
modelInd = irf_ask('Which shock model: [%]>','modelInd',1);

%% approved data indices

% print progress
fprintf('\n Setting parameters and limiting selection...\n')

% Solar wind energy array in [eV]
EswV = .5*u.mp*(sum((VuV*1e3).^2,2)/u.e);

% accepted intervals
dind = find((EfpiMaxV)./EswV>EmaxFac);

% all
%dind = 1:N;

% redefine N
N = numel(find(dind));

fprintf(['Number of crossings within parameters: ',num2str(N),' out of ',num2str(length(TV)),'\n'])

%% Clean data arrays
% single values
TV = TV(dind);
EswV = EswV(dind);
dTV = dTV(dind);
TuV = TuV(dind);
TdV = TdV(dind);
dTuV = dTuV(dind);
dTdV = dTdV(dind);
VuV = VuV(dind,:);
VuLV = VuLV(dind,:);
VdV = VdV(dind,:);
BuV = BuV(dind,:);
BuLV = BuLV(dind,:);
BdV = BdV(dind,:);
NuV = NuV(dind);
NuLV = NuLV(dind);
NdV = NdV(dind);
thBrV = thBrV(dind);
betaiV = betaiV(dind);
EV = EV(dind,:);
dEV = dEV(dind,:);
fdV = fdV(dind,:);
fuV = fuV(dind,:);
EfpiMaxV = EfpiMaxV(dind);
RV = RV(dind,:);
hasEISV = hasEISV(dind);
% convert to boolean
hasEISV = (hasEISV==1);
dstV = dstV(dind);
kpV = kpV(dind);
ssnV = ssnV(dind);
s107V = s107V(dind);
aeV = aeV(dind);
lineNumV = lineNumV(dind);
swIdV = swIdV(dind);
imfIdV = imfIdV(dind);

for jj = 1:length(shModel)
    nvecV.(shModel{jj}) = nvecV.(shModel{jj})(dind,:);
    MaV.(shModel{jj}) = MaV.(shModel{jj})(dind);
    MfV.(shModel{jj}) = MfV.(shModel{jj})(dind);
    thBnV.(shModel{jj}) = thBnV.(shModel{jj})(dind);
    thVnV.(shModel{jj}) = thVnV.(shModel{jj})(dind);
    sigV.(shModel{jj}) = sigV.(shModel{jj})(dind);
end

Nevents = numel(find(dind));

% angle between earth-sun line and sc position in xy plane
[alphaV,~] = cart2pol(RV(:,1),RV(:,2),RV(:,3));
alphaV = alphaV*180/pi; % degrees

% angle between earth-sun line and sc position
[phiV,~] = cart2pol(RV(:,1),sqrt(RV(:,2).^2+RV(:,3).^2));
phiV = phiV*180/pi; % degrees


%% Select shock model

fprintf(['Using bow shock model: ',shModel{modelInd},'\n'])
thBnV1 = thBnV.(shModel{modelInd});
thVnV1 = thVnV.(shModel{modelInd});
MaV1 = MaV.(shModel{modelInd});
MfV1 = MfV.(shModel{modelInd});
nvecV1 = nvecV.(shModel{modelInd});
sigV1 = sigV.(shModel{modelInd});


%% Ask if user wishes to calculate shock connenction time
fprintf('...done \n\n')
if optNum ~= 7 && optNum ~= 12
    calcTconn = irf_ask('Calculate shock age from model (slow)? (0:no 1:yes): [%]>','calcTconn',1);
    
    if calcTconn
        TcLim = irf_ask('Set upper limit of shock age in [s]: [%]>','TcLim',600);
        getShockAges;
    end
end
%% Calculate energy densities and define accelaration efficiency
[Ud0,Ud1,Ud2,Ud3] = get_energy_dens(EV,dEV,fdV,VuV,[3,5,10]);
[Uu0,Uu1,Uu2,Uu3] = get_energy_dens(EV,dEV,fuV,VuV,[3,5,10]);
accEffV1 = Ud3./Ud0;


%% calculate spectral slope of distribtion functions (limits are somewhat arbitrary)
specSlope = get_spectral_slope(EV,fdV,VuV,[5,15]);


