
% initialize structure
TcV = [];
nt = 2e2;
% for im = 1:length(shModel) % all models (not working now)
for im = modelInd % just the selected model (should be farris)
    TcV.(shModel{im}) = zeros(N,1);
    thBnHistoryV.(shModel{im}) = zeros(N,nt);
end


% ion cyclotron frequency
wcpV = zeros(N,1);

disp('Calculating shock connection time...')
fprintf('id = %4.0f/%4.0f\n',0,N) % display progress
for id = 1:N
    % display progress
    if mod(id,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],id,N); end 
    
    scd = [];
    scd.Bu = BuV(id,:);
    scd.Vu = VuV(id,:);
    scd.nu = NuV(id);
    scd.R = RV(id,:);
    
    
    %for im = 1:length(shModel) % all models (not working now)
    for im = modelInd % just the selected model (should be farris)
         bstemp = sh_bowshock_time(scd,'Nbs',500,'debug',0,'thMax',140,'shmodel',shModel{im},'nt',nt);
         bstemp.thBn(bstemp.thBn>90) = 180-bstemp.thBn(bstemp.thBn>90);
         TcV.(shModel{im})(id) = bstemp.t(end);
         TcValt.(shModel{im})(id) = bstemp.tend;
         % save all thBn vectors for debug reasons
         thBnHistoryV.(shModel{im})(id,:) = bstemp.thBn';
    end
    

    
%     Tl = bsspec.t(bsspec.thBn<60);
%     if ~isempty(Tl)
%         TcV(id) = bsspec.t(end)-Tl(1);
%     end
    
    
    wcpV(id) = norm(BuV(id,:))*1e-9*u.e/u.mp;
    
end



%%
% unly use the selected shock model
TcV1 = TcV.(shModel{modelInd});

% some special treatment
% only consider times below some value
idTc = find(TcV1<TcLim); % also gets rid of nans






% %% scatter plot
% 
% fig = figure;
% hca = axes(fig);
% 
% % plot events with Ma as color
% scatter(hca,phiV(idTc),TcV1(idTc).*wcpV(idTc),400,thBnV1(idTc),'.')
% hold(hca,'on')
% hca.XLim = [0,100];
% hca.YLim(1) = 0;
% 
% %sh_cmap(hca,cmap)
% colormap(jet)
% hca.Color = axcol;
% fig.Color = figcol;
% hca.XAxis.Color = textcol;
% hca.YAxis.Color = textcol;
% %grid(hca,'on')
% 
% hcb = colorbar(hca);
% hcb.Color = textcol;
% hca.CLim = [0,90];
% 
% hca.Box = 'on';
% 
% ylabel(hca,'$T_c$ [$\omega_{cp}^{-1}$]','interpreter','latex')
% xlabel(hca,'$\phi$ [$^{\circ}$]','interpreter','latex')
% ylabel(hcb,'$\theta_{Bn}$','Fontsize',17,'interpreter','latex')
% 
% hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
% hleg.BackgroundColor = hca.Color;
% 
% hca.LineWidth = 1.2;
% hca.FontSize = 17;
% hcb.LineWidth = 1.2;
% fig.InvertHardcopy = 'off';
% 




