function [] = sh_panel_span(ax,yl)
% SH_PANEL_SPAN Span panels in figure

N = numel(ax);

h = diff(yl)/N;

for k = 1:N
    ax(k).Position(2) = yl(1)+(N-k)*h;
    ax(k).Position(4) = h;
end

end



