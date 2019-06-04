function [hl] = plot_acceff_dep(ax,thBnV,thBinEdges,accEffV,XV,Xedges,colmat,varstr,unitstr,plotType)
%PLOT_ACCEFF_DEP Plots acceleration efficiency as a function of thBn
%
%   PLOT_ACCEFF_DEP(ax,thBnV,thBninEdges,accEffV,XV,Xedges,colmat,varstr,unitstr,plotType)
%       scatter plot of thBn and accEff with averages in bins defined by
%       thBinEdges. The averages are divided by some parameter given by XV
%       into classes defined by Xedges.
% 
%       Other input:
%           colmat  -   color matrix containing colors for classes
%                       col=colmat(1,:),...
%           varstr  -   string in LaTeX format for variable in XV
%           unitstr -   string in LaTeX format for unit of XV
%
%   hl = PLOT_ACCEFF_DEP(...) returns handle to legend in plot
%
% TODO: A lot

confval = 0.9;

dthBin = median(diff(thBinEdges));

idTh = discretize(thBnV,thBinEdges);
idX = discretize(XV,Xedges);

accEffAvg2D = zeros(length(thBinEdges)-1,length(Xedges)-1);
accEffStd2D = zeros(length(thBinEdges)-1,length(Xedges)-1);
accEffErrUp2D = zeros(length(thBinEdges)-1,length(Xedges)-1);
accEffErrDown2D = zeros(length(thBinEdges)-1,length(Xedges)-1);

for ii = 1:length(thBinEdges)-1
    for jj = 1:length(Xedges)-1
        
        % mean and std
        mu =  nanmean(accEffV(idTh==ii & idX==jj));
        sig = nanstd(accEffV(idTh==ii & idX==jj));
        
        % after a lot(!) of math
        accEffErrUp2D(ii,jj) = norminv(confval*normcdf(mu/sig)-normcdf(mu/sig)+1)*sig;
        accEffErrDown2D(ii,jj) = norminv(1-confval*normcdf(mu/sig))*sig;
        
        accEffAvg2D(ii,jj) = mu;
        accEffStd2D(ii,jj) = sig;
       
    end
end


thoffs = [-1,0,1,2];


%% Plot
hold(ax,'on')
switch lower(plotType)
    
    case 'line'
        %% plot line
        
        
        
        
        for jj = 1:size(accEffAvg2D,2)
            scatter(ax,thBnV(idX==jj),accEffV(idX==jj)*100,40,colmat(jj,:),'o','MarkerFaceColor','flat');
        end
        
        
        for jj = 1:size(accEffAvg2D,2)
            errorbar(ax,thBinEdges(1:end-1)+dthBin/2+thoffs(jj),accEffAvg2D(:,jj)*100,accEffErrDown2D(:,jj)*100,accEffErrUp2D(:,jj)*100,'linewidth',3,'color',colmat(jj,:))
        end
        
        
        legStr = cell(1,size(accEffAvg2D,2));
        for jj = 1:size(accEffAvg2D,2)
            if jj == 1
                legStr{jj} = ['$',varstr,'{<}',num2str(Xedges(jj+1)),'$\,',unitstr,', \hspace{.5cm} $N = ',num2str(numel(find(idX==jj))),'$'];
            elseif jj == size(accEffAvg2D,2)
                %unitstr = unitstrtemp;
                legStr{jj} = ['$',varstr,'{>}',num2str(Xedges(jj)),'$\,',unitstr,', \hspace{.5cm} $N = ',num2str(numel(find(idX==jj))),'$'];
            else
                % hack
                unitstrtemp = unitstr;
                unitstrtemp(ismember(unitstrtemp,'$')) = '';
                legStr{jj} = ['$',num2str(Xedges(jj)),'\,',unitstrtemp,'{<}',varstr,'{<}',num2str(Xedges(jj+1)),'\,',unitstrtemp,'$, \hspace{.5cm} $N = ',num2str(numel(find(idX==jj))),'$'];
                
            end
        end
        
        
        hl = legend(ax,legStr);
        hl.Interpreter = 'latex';
        hl.FontSize = 15;
        
        ylabel(ax,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')

        
    case 'hist'
        
        surf(thBinEdges,Xedges,zeros(length(Xedges),length(thBinEdges)),accEffAvg2D'*100)
        ylabel(ax,['$',varstr,'$ [',unitstr,']'],'Fontsize',15,'interpreter','latex')
        view(2)
        scatter(thBnV,XV,200,'k.')
        hcb = colorbar(ax);
        ylabel(hcb,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
end

ax.Box = 'on';

ax.XLim = [0,90];
%hca.YLim = [0,15];
ax.YLim(1) = 0;

xlabel(ax,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

%irf_legend(hca,[num2str(beta*100),'\% significance'],[0.98,0.92],'Fontsize',15,'interpreter','latex','color',textcol);
%hl.BackgroundColor = hca.Color;

ax.LineWidth = 1.2;
ax.FontSize = 14;


end

