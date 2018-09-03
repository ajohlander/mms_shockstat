function [out] = sh_cmap(varargin)
% SH_CMAP Custom colormap.
%
%   SH_CMAP(cmp_name) Applies the given colormap to current axes.
%
%   SH_CMAP(AX,...) Applies the given colormap to axes AX.
%
%   cmp = SH_CMAP(AX,...) Returns the colormap.
%
%   Colormap Names:
%       'standard'  -   similar to irf_colormap
%       'plopp'      -   orange and blue
%       'bluered'   -   like irf but white in middle
%       'speedway'  -   white-yellow-red-blue
%
%   Colors:
%       'red'
%       'green'
%       'blue'
%       'yellow'
%       'orange'
%       'lime'
%       'gray'
%
%   See also: IRF_COLORMAP

%% Input
[ax,args,nargs] = axescheck(varargin{:});

% check which axis to apply
if isempty(ax)
    axes(gca);
end

N = 64; % number of points
if(nargs == 0)
    cMapMode = 'standard';
else
    cMapMode = args{1};
end

didFindColormap = 1;

%% Set which colors
switch lower(cMapMode) % Special colormaps
    case 'standard'
        c = [255,255,255;...
            043,255,255;...
            000,255,000;...
            255,255,000;...
            255,255,000;...
            255,000,000;...
            000,000,255]/255;
    case 'irf'
        c = [255,255,255;...
            043,255,255;...
            000,255,000;...
            255,255,000;...
            255,000,000;...
            000,000,255]/255;
    case 'plopp'
        c = [1,1,1;...
            1,0.7,0;...
            0,0,0.8;...
            0,0,0];
    case 'bluered'
        c = [0,0,1;...
            1,1,1;...
            1,0,0];
    case 'redblue'
        c = [1,0,0;...
            1,1,1;...
            0,0,1];
    case 'speedway'
        c = [255,255,255;...
            255,255,000;...
            255,000,000;...
            000,000,200]/255;
    case 'strangeways'
        c = [55,137,187;...
            106,193,165;...
            172,220,166;...
            230,244,157;...
            255,254,194;...
            253,223,144;...
            251,173,104;...
            242,109,074;...
            211,064,082]/255;
    otherwise % Single color
        
        switch lower(cMapMode)
            case 'red'
                cMiddle = [1,0,0];
            case 'green'
                cMiddle = [0,1,0];
            case 'blue'
                cMiddle = [0,0,1];
            case 'yellow'
                cMiddle = [1,1,0];
            case 'orange'
                cMiddle = [1,0.7,0];
            case 'lime'
                cMiddle = [0.5,1,0];
            case 'gray'
                cMiddle = [0.5,0.5,0.5];
            case 'grey'
                cMiddle = [0.5,0.5,0.5];
            otherwise
                didFindColormap = 0;
        end
        if didFindColormap
            c = [1,1,1;cMiddle;0,0,0];
        end
end

if ~didFindColormap
    % try Matlabs own colormaps
    try
        eval(['c = ',cMapMode,';'])
    catch
        error('Unknown color')
    end
end

%% Create colormap
cN = size(c,1);
x = linspace(1,N,cN);
cmap = zeros(N,3);

for i = 1:N
    cmap(i,:) = interp1(x,c,i);
end

for ii = 1:length(ax)
    colormap(ax(ii),cmap);
end

if nargout > 0
    out = cmap;
end

end

