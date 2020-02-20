function [] = sh_panel_span(ax,yl)
% SH_PANEL_SPAN Span panels in figure

N = numel(ax);

xl = [];
if length(yl)==4
    xl = yl([1,3]);
    yl = yl([2,4]);
end

h = diff(yl)/N;
l = diff(xl);

for k = 1:N
    ax(k).Position(2) = yl(1)+(N-k)*h;
    ax(k).Position(4) = h;
end

if ~isempty(xl)
    for k = 1:N
        ax(k).Position(1) = xl(1);
        ax(k).Position(3) = l;
    end
end

end



