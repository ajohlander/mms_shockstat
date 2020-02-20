function varargout = sh_figure(subNum,figSize,varargin)
% SH_FIGURE Quick way to call irf_plot.
%
%   SH_FIGURE(subNum) initiates figure with number of panels subNum.
%   Size of figure is set to 10x7 cm.
%
%   SH_FIGURE(subNum,figSize) declare size of figure with a vector
%   figSize = [width, height] in centimeters.
%   h = SH_FIGURE(...) returns vector of axes for the panels h
%
%   [h,f] = SH_FIGURE(...) also returns figure f.
%
%   See also: IRF_PLOT


%% Input
if(nargin == 1)
    figSize = [10,7];
elseif(nargin == 0)
    figSize = [10,7];
    subNum = 1;
end


%% Check for flags

% Set default values
regularPlot = 0;

if nargin == 3 && strcmpi(varargin{1},'noirf')
    regularPlot = 1;
end
    

%% Initiate figures
% Initiate figure
irf.log('w','Initiating new figure')

if regularPlot
    f = figure;
    f.Color = 'w';
    h = gobjects(1,subNum);
    for ii = 1:subNum
        h(ii) = subplot(subNum,1,ii);
        h(ii).Box = 'on';
    end
    linkaxes(h,'x')
    sh_panel_span(h,[.1,.95])
    f.UserData.subplot_handles = h;
else
    h = irf_plot(subNum,'newfigure');
    f = h.Parent;
end

% Set parameters
f.PaperUnits = 'centimeters';
xSize = figSize(1); ySize = figSize(2);
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
f.PaperPosition = [xLeft yTop xSize ySize];
f.Position = [10 10 xSize*50 ySize*50];
f.PaperPositionMode = 'auto';

for i = 1:subNum
    h(i).UserData.XLabel.String = h(i).XLabel.String; % Stores XLabel for restoration
    h(i).Layer = 'top'; % Axes visible over patches (and surfaces?).
    if(i ~= subNum)
        h(i).UserData.XLabel.Visible = 'off';
        h(i).XTickLabel = '';
        h(i).XLabel.String = '';
    else
        h(i).UserData.XLabel.Visible = 'on';
    end
end



%% Output
if nargout == 1
    varargout{1} = h;
elseif nargout == 2
    varargout{1} = h;
    varargout{2} = f;
end


