% Compare shock models


%% scatter plots combined with chosen model


fig = figure;


axCount = 1;
axInd = 1:length(shModel); axInd(modelInd) = [];

for jj = axInd
    
    hca = subplot(2,3,axCount);
    
    hold(hca,'on')
    
    scatter(hca,thBnV1,thBnV.(shModel{jj})-thBnV1,300,alphaV,'.')
    
    plot(hca,[0,90],[0,0],'k--','linewidth',1.3)
    
    xlabel(hca,['$\theta_{',shModel{modelInd},'}$ [$^{\circ}$]'],'Fontsize',15,'interpreter','latex')
    ylabel(hca,['$\theta_{',shModel{jj},'}-\theta_{',shModel{modelInd},'}$ [$^{\circ}$]'],'Fontsize',15,'interpreter','latex')
    axCount = axCount+1;
    
    axis(hca,'equal')
    box(hca,'on')
    hca.XLim = [0,90]; hca.YLim = [-45,45];
    sh_cmap(hca,cmap)
    hca.CLim = [0,90];
    
    hca.Color = axcol;
    hca.XAxis.Color = textcol;
    hca.YAxis.Color = textcol;    
    hca.LineWidth = 1.2;
end

axPos = hca.Position;
hcb = sh_cbar(hca);
hca.Position = axPos;
ylabel(hcb,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

hcb.Position(1) = sum(hca.Position([1,3]))+.1;
fig.Color = figcol;
hcb.Color = textcol;

hcb.LineWidth = 1.2;




%% compression factors and angle to sun-earth line

fig = figure;
hca = axes(fig);
hold(hca,'on')

for jj = 1:length(shModel)
    plot(hca,alphaV,sigV.(shModel{jj}),'.','markersize',15)
end
plot(hca,[0,90],[1,1],'k--','linewidth',1.3)

legStr = cell(1,length(shModel));
for jj = 1:length(shModel)
    legStr{jj} = [shModel{jj},':   $',...
        num2str(round(mean(sigV.(shModel{jj})),2)),'\pm'...
        num2str(round(std(sigV.(shModel{jj})),2)),'$'];
end

hleg = legend(legStr);
hleg.FontSize = 15;
hleg.Interpreter = 'latex';

hca.XLim = [0,max(alphaV)]; hca.YLim = [.5,1.5];

xlabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
ylabel(hca,'$\sigma$','Fontsize',15,'interpreter','latex')


fig.Color = figcol;
hca.Color = axcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
hca.LineWidth = 1.2;
box(hca,'on')



%% "standard deviation" of normal angle to average normal vs alpha

% do math

nStdAng = zeros(length(nvecV.(shModel{1})),1);
nAvg = zeros(length(nvecV.(shModel{1})),3);

for ii = 1:length(nvecV.(shModel{1}))
    nTemp = zeros(length(shModel),3);
    for jj = 1:length(shModel)
        nTemp(jj,:) = nvecV.(shModel{jj})(ii,:);
    end
    nAvg(ii,:) = mean(nTemp);
    nAvg(ii,:) = nAvg(ii,:)/norm(nAvg(ii,:));
    
    nAngTemp = zeros(length(shModel),1);
    for jj = 1:length(shModel)
        % not really std
        nAngTemp(jj) = acosd(dot(nvecV.(shModel{jj})(ii,:),nAvg(ii,:)));
    end
    nStdAng(ii,:) = mean(nAngTemp);
    
end


fig = figure;
hca = axes(fig);
hold(hca,'on')

plot(alphaV,nStdAng,'.','color',col1,'markersize',25)

hca.XLim = [0,max(alphaV)]; hca.YLim(1) = 0;

ylabel(hca,'$\langle \angle \mathbf{n}_i, \langle \mathbf{n}\rangle \rangle $ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

fig.Color = figcol;
hca.Color = axcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
hca.LineWidth = 1.2;
box(hca,'on')





%% maximum deviation of normal angle to average normal vs alpha

% do math

nMaxAng = zeros(length(nvecV.(shModel{1})),1);

for ii = 1:length(nvecV.(shModel{1}))
    
    nAngTemp = zeros(length(shModel),1);
    for jj = 1:length(shModel)
        for kk = 1:length(shModel)
            angTemp = acosd(dot(nvecV.(shModel{jj})(ii,:),nvecV.(shModel{kk})(ii,:)));
            if angTemp > nMaxAng(ii,:)
                nMaxAng(ii,:) = angTemp;
            end
        end
    end

    
end


fig = figure;
hca = axes(fig);
hold(hca,'on')

plot(alphaV,nMaxAng,'.','color',col1,'markersize',25)

ylabel(hca,'$\max{ (\angle \mathbf{n}_i, \mathbf{n}_j)} $ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

hca.XLim = [0,max(alphaV)]; hca.YLim(1) = 0;

fig.Color = figcol;
hca.Color = axcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
hca.LineWidth = 1.2;
box(hca,'on')





