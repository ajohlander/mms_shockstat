%% read time intervals

rePlotAll = irf_ask('Redo all plots? (0:no, 1:yes) [%]>','rePlotAll',0);

filePath = irf_ask('Choose path for plots [%]>','filePath','../mms_shock_plots/mmsVomni/');

edgeMinutes = irf_ask('Add how many minutes on each side of plot? [%]>','edgeMinutes',5);

u = irf_units;

%% plot position parameters
axl = 0.125; % axis left
axu = 0.06;
axw = .82;
axh = .91;
cbl = axl+axw+.01;
cbw = .03;


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
    tint = tint+edgeMinutes*60*[-1,1];
    
    
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
    B = mms.get_data('B_gse_fgm_srvy_l2',tint,ic);
    
    if isempty(B)
        disp(['no FGM data for ',tstr,', skipping...'])
        continue;
    end
    
    ni = mms.get_data('Ni_fpi_fast_l2',tint,ic);
    ne = mms.get_data('Ne_fpi_fast_l2',tint,ic);
    Vi = mms.get_data('Vi_gse_fpi_fast_l2',tint,ic);
    
    % check if all data is there
    if isempty(B) || isempty(ni) || isempty(ne) || isempty(Vi)
        disp(['no data for ',tstr,', skipping...'])
        continue;
    end
    
    
    %% read omni data
    
    ff= irf_get_data(tint,'bx,by,bz,vx,vy,vz,n','omni_min');
    
    if isempty(ff)
        irf.log('c','failed to read omni data, skipping...')
        break;
    else
        omniTime = irf_time(ff(:,1),'epoch>epochTT');
        Bomni = irf.ts_vec_xyz(omniTime,ff(:,2:4));
        Vomni = irf.ts_vec_xyz(omniTime,ff(:,5:7)+[0,29.8,0]);
        Nomni = irf.ts_scalar(omniTime,ff(:,8));
    end

    
    %% plot
    h = sh_figure(8,[12,16]);
    
    % color order
    col = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250];
    
    % babs
    hca = irf_panel(h,'babs');
    hold(hca,'on')
    irf_plot(hca,Bomni.abs,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,B.abs,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B$ [nT]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    hca.YLim(1) = 0;
    hleg = irf_legend(hca,['MMS',num2str(ic)],[0.02,0.95],'Fontsize',15,'interpreter','latex');
    hleg.BackgroundColor = 'w';
    hca.ColorOrder = col;
    % just make regular legend
    irf_legend(hca,{'FPI-DIS','OMNI','FPI-DES'},[0.02,1.04],'Fontsize',15,'interpreter','latex')
    
    % bx
    hca = irf_panel(h,'bx');
    hold(hca,'on')
    irf_plot(hca,Bomni.x,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,B.x,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_x$ [nT]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    % by
    hca = irf_panel(h,'by');
    hold(hca,'on')
    irf_plot(hca,Bomni.y,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,B.y,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_y$ [nT]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    % bz
    hca = irf_panel(h,'bz');
    hold(hca,'on')
    irf_plot(hca,Bomni.z,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,B.z,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$B_z$ [nT]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    
    % vx
    hca = irf_panel(h,'vx');
    hold(hca,'on')
    irf_plot(hca,Vomni.x,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,Vi.x,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_x$ [km/s]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    % vy
    hca = irf_panel(h,'vy');
    hold(hca,'on')
    irf_plot(hca,Vomni.y,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,Vi.y,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_y$ [km/s]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    % vy
    hca = irf_panel(h,'vz');
    hold(hca,'on')
    irf_plot(hca,Vomni.z,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,Vi.z,'color',col(1,:),'linewidth',1.5);
    ylabel(hca,'$V_z$ [km/s]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    % n
    hca = irf_panel(h,'n');
    hold(hca,'on')
    irf_plot(hca,Nomni,'-x','color',col(2,:),'linewidth',3);
    irf_plot(hca,ni,'color',col(1,:),'linewidth',1.5);
    irf_plot(hca,ne,'color',col(3,:),'linewidth',1.5);
    ylabel(hca,'$N$ [cm$^{-3}$]','fontsize',15,'interpreter','latex')
    hca.YLimMode ='auto';
    
    
    irf_zoom(h,'x',tint)
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
        
        irf_pl_mark(h(jj),tint+edgeMinutes*60*[1,-1])
    end
    
    hcb1.LineWidth = 1.3; hcb2.LineWidth = 1.3;
    hca = irf_panel(h,'fi');
    hcb1.Position([2,4]) = hca.Position([2,4]);
    hcb1.Position([1,3]) = [cbl,cbw];
    hca = irf_panel(h,'redi');
    hcb2.Position([2,4]) = hca.Position([2,4]);
    hcb2.Position([1,3]) = [cbl,cbw];
    
    
    %% save figure
    % give it the time stamp of the start time
    irf_print_fig(h(1),[filePath,fileName(1:end-4)],'png')
    
    %% close figure
    close(h(1).Parent)
    
end


