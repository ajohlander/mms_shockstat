%% read time intervals

rePlotAll = irf_ask('Redo all plots? (0:no, 1:yes) [%]>','rePlotAll',0);

filePath = irf_ask('Choose path for plots [%]>','filePath','../mms_shock_plots/mmsVomni/');

edgeMinutes = irf_ask('Add how many minutes on each side of plot? [%]>','edgeMinutes',5);

u = irf_units;

%% plot position parameters
axl = 0.125; % axis left
axu = 0.06;
axw = .82;
axh = .89;
cbl = axl+axw+.01;
cbw = .03;


%% Set limits in y direction in plots
% factor times norm(Bu)
fB = 2;
% plus/minus Bu(i)
dB = 5;
% plus/minus Vu(i)
dV = 100;
% factor times Nu
fN = 2;

%% Make plot
tline = 1;
lineNum = 0;

fid = fopen(listFileName);

% loop through skipped lines
for ii = 1:startLine-1
    lineNum = lineNum+1;
    tline = fgets(fid);
end

% super trooper looper
while tline ~= -1
    lineNum = lineNum+1;
    disp(['Current line number: ',num2str(lineNum)])
    
    %% end loop if requested
    if lineNum >= stopLine; disp('Reached stop line, exiting...'); break; end
    
    %% read line from file
    tline = fgets(fid);
    
    if tline(1) == -1 || ~strcmp(tline(1),'2')
        disp('done!')
        break;
    end
    tintStr = [tline(1:10),'T',tline(12:19),'/',tline(23:32),'T',tline(34:41)];
    
    % all other strings if needed
    restStr = tline(44:end);
    
    % find commas
    idcomma = strfind(restStr,',');
    % find start paranthesis
    idpar = strfind(restStr,'(');
    
    fomStr = restStr(1:idcomma(1)-1);
    sitlStr = restStr(idcomma(1)+2:idpar(1)-1);
    
    %%
    descStr = restStr(idcomma(2)+2:end);
    if length(descStr)>50
        descStr = descStr(1:50);
    end
    
    %% get time interval and check if already done
    tint = irf.tint(tintStr);
    
    % add some time 
    tint2 = tint+edgeMinutes*60*[-1,1];
    
    
    tstr = tint(1).toUtc;
    tstr([5,8,14,17]) = '';
    tstr = tstr(1:15);
    
    listing = dir(filePath);
    fnCell = {listing.name};
    fileName = ['sh_mmsVomni_',tstr,'.png'];
    
    if ismember(fileName,fnCell) && ~rePlotAll
        disp(['Skipping ',tstr])
        continue;
    end
    
    
    
    %% get data 
    B = mms.get_data('B_gse_fgm_srvy_l2',tint2,ic);
    
    if isempty(B)
        disp(['no FGM data for ',tstr,', skipping...'])
        continue;
    end
    
    ni = mms.get_data('Ni_fpi_fast_l2',tint2,ic);
    ne = mms.get_data('Ne_fpi_fast_l2',tint2,ic);
    Vi = mms.get_data('Vi_gse_fpi_fast_l2',tint2,ic);
    
    % check if all data is there
    if isempty(B) || isempty(ni) || isempty(ne) || isempty(Vi)
        disp(['no data for ',tstr,', skipping...'])
        continue;
    end
    
    
    %% read omni data
    
    ff= irf_get_data(tint2,'bx,by,bz,vx,vy,vz,n','omni_min');
    
    if isempty(ff)
        irf.log('c','failed to read omni data, skipping...')
        break;
    else
        omniTime = irf_time(ff(:,1),'epoch>epochTT');
        Bomni = irf.ts_vec_xyz(omniTime,ff(:,2:4));
        % correct velocity for abberation 
        Vomni = irf.ts_vec_xyz(omniTime,ff(:,5:7)+[0,29.8,0]);
        Nomni = irf.ts_scalar(omniTime,ff(:,8));
    end
    
    %for good measures, remove old ff
    clear ff

    
    %% get average upstream omni values
    
    Bu = nanmean(Bomni.tlim(tint).data,1);
    if isnan(Bu); Bu = nan(1,3); end
    Vu = nanmean(Vomni.tlim(tint).data,1);
    if isnan(Vu); Vu = nan(1,3); end
    nu = nanmean(Nomni.tlim(tint).data,1);
    
    %% plot
    h = sh_figure(8,[12,16]);
    
    % color order
    col = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250];
    
    % babs
    hca = irf_panel(h,'babs');
    hold(hca,'on')
    irf_plot(hca,Bomni.abs,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,norm(Bu)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,B.abs,'color',col(1,:),'linewidth',1.5);
    % for legend
    plot(hca,nan,'color',col(3,:),'linewidth',1.5)
    ylabel(hca,'$B$ [nT]','fontsize',15,'interpreter','latex')
    %hca.YLimMode ='auto';
    if ~isnan(Bu(1)); hca.YLim = [0,fB*norm(Bu)]; end
    hleg = irf_legend(hca,['MMS',num2str(ic)],[0.02,0.95],'Fontsize',15,'interpreter','latex');
    hleg.BackgroundColor = 'w';
    hca.ColorOrder = col;
    % make legend
    hleg = legend(h(1),'OMNI','Upstream','FGM/DIS','DES','location','northoutside');
    hleg.Orientation = 'horizontal';
    hleg.Position(2) = 0.947;
    hleg.Interpreter = 'latex';
    hleg.FontSize = 15;
    hleg.Box = 'off';
    
    
    % bx
    hca = irf_panel(h,'bx');
    hold(hca,'on')
    irf_plot(hca,Bomni.x,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Bu(1)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,B.x,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_x$ [nT]','fontsize',15,'interpreter','latex')
    if ~isnan(Bu(1)); hca.YLim = Bu(1)+[-1,1]*dB; end
    
    % by
    hca = irf_panel(h,'by');
    hold(hca,'on')
    irf_plot(hca,Bomni.y,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Bu(2)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,B.y,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_y$ [nT]','fontsize',15,'interpreter','latex')
    if ~isnan(Bu(2)); hca.YLim = Bu(2)+[-1,1]*dB; end
    
    % bz
    hca = irf_panel(h,'bz');
    hold(hca,'on')
    irf_plot(hca,Bomni.z,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Bu(3)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,B.z,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_z$ [nT]','fontsize',15,'interpreter','latex')
    if ~isnan(Bu(3)); hca.YLim = Bu(3)+[-1,1]*dB; end
    
    % vx
    hca = irf_panel(h,'vx');
    hold(hca,'on')
    irf_plot(hca,Vomni.x,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Vu(1)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,Vi.x,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_x$ [km/s]','fontsize',15,'interpreter','latex')
    if ~isnan(Vu(1)); hca.YLim = Vu(1)+[-1,1]*dV; end
    
    % vy
    hca = irf_panel(h,'vy');
    hold(hca,'on')
    irf_plot(hca,Vomni.y,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Vu(2)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,Vi.y,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_y$ [km/s]','fontsize',15,'interpreter','latex')
    if ~isnan(Vu(2)); hca.YLim = Vu(2)+[-1,1]*dV; end
    
    % vy
    hca = irf_panel(h,'vz');
    hold(hca,'on')
    irf_plot(hca,Vomni.z,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,Vu(3)*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,Vi.z,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_z$ [km/s]','fontsize',15,'interpreter','latex')
    if ~isnan(Vu(3)); hca.YLim = Vu(3)+[-1,1]*dV; end
    
    % n
    hca = irf_panel(h,'n');
    hold(hca,'on')
    irf_plot(hca,Nomni,'o','color',col(2,:),'linewidth',1.8,'markersize',6);
    irf_plot(hca,[tint2.epochUnix,nu*[1;1]],'.-','color',col(2,:),'linewidth',1.5);
    irf_plot(hca,ni,'color',col(1,:),'linewidth',1.5);
    irf_plot(hca,ne,'color',col(3,:),'linewidth',1.5);
    ylabel(hca,'$N$ [cm$^{-3}$]','fontsize',15,'interpreter','latex')
    if ~isnan(nu); hca.YLim = [0,fN*norm(nu)]; end
    
    irf_zoom(h,'x',tint2)
    irf_plot_axis_align(h)
    sh_panel_span(h,[axu,axu+axh])
    pause(0.01)
    for jj = 1:length(h)
        irf_zoom(h(jj),'y',h(jj).YLim)
        h(jj).Position(1) = axl;
        h(jj).Position(3) = axw;
        
        h(jj).Layer = 'top';
        h(jj).FontSize = 15;
        h(jj).LineWidth = 1.3;
        h(jj).YLabel.Position(1) = -0.08;
        
        irf_pl_mark(h(jj),tint)
    end
    
    hcb1.LineWidth = 1.3; hcb2.LineWidth = 1.3;
    hca = irf_panel(h,'fi');
    hcb1.Position([2,4]) = hca.Position([2,4]);
    hcb1.Position([1,3]) = [cbl,cbw];
    hca = irf_panel(h,'redi');
    hcb2.Position([2,4]) = hca.Position([2,4]);
    hcb2.Position([1,3]) = [cbl,cbw];
    
    %% make title
    % title with file name and line number
    irf_legend(h(1),tstr,[0,1.3],'Fontsize',15,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['Line number: ',num2str(lineNum)],[1,1.3],'Fontsize',15,'interpreter','latex','horizontalalignment','right')
    
    
    %% save figure
    % give it the time stamp of the start time
    irf_print_fig(h(1),[filePath,fileName(1:end-4)],'png')
    
    %% close figure
    close(h(1).Parent)
    
end


