function FDthr = tmfc_plot_FD(FD,options,SPM_paths,subject_paths,struct_paths,funct_paths)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Opens a GUI with FD plots. Allows the user to select the FD threshold
% for spike regression.
%
% -------------------------------------------------------------------------
% FORMAT: FDthr = tmfc_plot_FD(FD)
% Allow to save group FD statistics only.
%
% FORMAT: FDthr = tmfc_plot_FD(FD,options,SPM_paths,subject_paths,struct_paths,funct_paths)
% Allow to save group FD statistics anf TMFC denoise settings.
% -------------------------------------------------------------------------
%
% =========================================================================
%
% Copyright (C) 2025 Ruslan Masharipov
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% Contact email: masharipov@ihb.spb.ru

% Display spike regression note
if nargin > 1
    if options.spikereg == 1
        spikenote = 'on';
    else
        spikenote = 'off';
    end
else
    spikenote = 'off';
end

% FD input only
if nargin == 1
    options = []; SPM_paths = []; subject_paths = []; struct_paths = []; funct_paths = [];
end

% Create GUI figure
FD_MW = figure('Name','Framewise displacement','NumberTitle','off','Units','normalized','Position',[0.20 0.1 0.55 0.80],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@FD_MW_exit);
FD_MW_txt1 = uicontrol(FD_MW,'Style','text','String','Select subject:','Units','normalized','Position',[0.075 0.94 0.85 0.038],'fontunits','normalized','FontSize',0.55,'HorizontalAlignment','Left','backgroundcolor','w');
FD_MW_LB1 = uicontrol(FD_MW,'Style','listbox','String',[],'Max',1,'Value',1,'Units','normalized','Position',[0.075 0.76 0.85 0.180],'FontUnits','points','FontSize',11.5,'callback',@update_plot);
movegui(FD_MW,'center');

% Create axes for plot
ax_frame = axes(FD_MW,'Units','normalized','Position',[0.075 0.43 .85 .3],'color',[0.75 0.75 0.75]);
box on; xlabel(ax_frame,'Scans','FontSize',9); ylabel(ax_frame,'FD, [mm]','FontSize',9); S = plot([0]);

% GUI elements
note_str = {'NOTE: Selected FD threshold will be used to create spike regressors (SpikeReg). For each flagged time point, a unit impulse function is included in general linear model. The number of spike regressors is equal to the number of flagged scans.'};
FD_MW_panel_1 = uipanel(FD_MW,'Units','normalized','Position',[0.075 0.15 0.48 0.22],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
FD_MW_panel_2 = uipanel(FD_MW,'Units','normalized','Position',[0.57 0.15 0.356 0.22],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
FD_MW_info_box = uicontrol(FD_MW,'Style','text','String',[],'Units','normalized','Position',[0.09 0.160 0.45 0.2],'fontunits','normalized','FontSize',0.09,'HorizontalAlignment','left','backgroundcolor','w');
FD_MW_threshold = uicontrol(FD_MW,'Style','pushbutton','String','FD threshold [mm]:','Units','normalized','Position',[0.66 0.27 0.18 0.050],'FontUnits','normalized','FontSize',0.34,'callback',@get_FDthr);
FD_MW_thres_edit = uicontrol(FD_MW,'Style','Edit','String','0.5','Units','normalized','Position',[0.66 0.20 0.18 0.050],'FontUnits','normalized','FontSize',0.4,'backgroundcolor','w');
FD_MW_note_txt = uicontrol(FD_MW,'Style','text','String',note_str,'Units','normalized','Position',[0.075 0.06 0.85 0.078],'fontunits','normalized','FontSize',0.22,'HorizontalAlignment','left','backgroundcolor','w','Visible',spikenote);
FD_MW_OK = uicontrol(FD_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.30 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@confirm_FDthr);
FD_MW_SAVE = uicontrol(FD_MW,'Style','pushbutton','String','Save','Units','normalized','Position',[0.55 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@save_data);

flagged = []; N_25prc = 0; N_50prc = 0; N_75prc = 0; N_95prc = 0;
mean_flagged = 0; sd_flagged = 0; max_flagged = 0; min_flagged = 0;

% Default FD threshold
%--------------------------------------------------------------------------
FDthr = 0.5;
plot_data(FDthr,1);

% Update FD plot
%--------------------------------------------------------------------------
function update_plot(~,~)
    selected_subject = get(FD_MW_LB1,'Value');
    plot_data(FDthr,selected_subject);
end

% Close GUI
%--------------------------------------------------------------------------
function FD_MW_exit(~,~)
    FDthr = 0.5;
    if spikenote == 1
        warning('The default FD threshold of 0.5 mm will be used for spike regression.');
    end
    delete(FD_MW);
end

% Change FDthr
%--------------------------------------------------------------------------
function get_FDthr(~,~)    
    temp_prefix = str2double(get(FD_MW_thres_edit,'String'));
    if isnan(temp_prefix)
        warning('Please enter a non-negative value for the FD threshold.');
    elseif temp_prefix < 0
        warning('Please enter a non-negative value for the FD threshold.');
    else
        FDthr = temp_prefix;
        selected_subject = get(FD_MW_LB1,'Value');
        plot_data(FDthr,selected_subject);
    end
end

% Export selected FDthr
%--------------------------------------------------------------------------
function confirm_FDthr(~,~)
    fprintf('Selected FD threshold is: %.3f mm.\n',FDthr);
    delete(FD_MW);
end

% Generate FD plot
%--------------------------------------------------------------------------
function plot_data(FD_threshold,iSub)
              
    delete([ax_frame,S]);
    
    % Create plot outline
    ax_frame = axes('Units','normalized','Position',[0.075 0.43 .85 .3],'color',[0.75 0.75 0.75]);
    box on; xlabel(ax_frame,'Scans','FontSize',9); ylabel(ax_frame,'FD, [mm]','FontSize',9);
    
    FD_ts = []; sess_sum = 0; sess = 0;
    for jSess = 1:length(FD(iSub).Sess)
        FD_ts = [FD_ts; FD(iSub).Sess(jSess).FD_ts];
        sess_sum = sess_sum + length(FD(iSub).Sess(jSess).FD_ts) + 1;
        sess = [sess; sess_sum];
    end

    % Plot FD time-series
    S = plot(FD_ts,'color',[0 0.447 0.7410 0.8]); xlabel('Scans'); ylabel('FD, [mm]'); xlim tight; x = xlim; y = ylim; hold on;
    % Plot sessions
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            plot([sess(jSess) sess(jSess)],[y(1) y(2)],'-k'); 
        end
        text(sess(jSess)+10,y(2),{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold on;
    end
    % Plot FD threshold
    S = plot([x(1) x(2)],[FD_threshold FD_threshold],'--','color',[0.9 0.2 0.2 0.8],'linewidth',2); text(x(2)+5,FD_threshold,{'FDthr'}); 
    update_text();
end

% Update text info
%--------------------------------------------------------------------------
function update_text(~,~)

    for jSub = 1:length(FD) 
        scans = 0;
        for jSess = 1:length(FD(jSub).Sess)
            flagged(jSub).Sess(jSess) = sum(FD(jSub).Sess(jSess).FD_ts > FDthr);
            scans = scans + length(FD(jSub).Sess(jSess).FD_ts);
        end
        flagged(jSub).total = sum(flagged(jSub).Sess);
        flagged(jSub).total_prc = 100*flagged(jSub).total/scans;
        clear scans
    end

    N_25prc = sum([flagged.total_prc]>25);
    N_50prc = sum([flagged.total_prc]>50);
    N_75prc = sum([flagged.total_prc]>75);
    N_95prc = sum([flagged.total_prc]>95);
    mean_flagged = round(mean([flagged.total]),1); mean_flagged_prc = round(mean([flagged.total_prc]),1);
    sd_flagged = round(std([flagged.total]),1); sd_flagged_prc = round(std([flagged.total_prc]),1);
    max_flagged = max([flagged.total]); max_flagged_prc = max([flagged.total_prc]);
    min_flagged = min([flagged.total]); min_flagged_prc = min([flagged.total_prc]);
    
    lb1_str = {};
    for jSub = 1:length(FD)
        lb1_str =  vertcat(lb1_str,strcat(num2str(FD(jSub).Subject),' :: (',num2str(round(flagged(jSub).total_prc,1),'%.1f'),'% above FDthr) :: (',num2str(flagged(jSub).total),' scans flagged) '));  
    end
    
    stat_str = {strcat(num2str(N_25prc),' subjects have >25% scans above FD threshold'),...
        strcat(num2str(N_50prc),' subjects have >50% scans above FD threshold'),...
        strcat(num2str(N_75prc),' subjects have >75% scans above FD threshold'),...
        strcat(num2str(N_95prc),' subjects have >95% scans above FD threshold'),...
        '------------------------------------------------------------------------------',...
        strcat('Mean number of flagged scans across subjects:',[' ' num2str(round(mean_flagged,1),'%.1f')],' (',num2str(round(mean_flagged_prc,1),'%.1f'),'%)'),...
        strcat('SD number of flagged scans across subjects:',  [' ' num2str(round(sd_flagged,1),'%.1f')],' (',num2str(round(sd_flagged_prc,1),'%.1f'),'%)'),...
        strcat('Max number of flagged scans across subjects:', [' ' num2str(max_flagged)],' (',num2str(round(max_flagged_prc,1),'%.1f'),'%)'),...
        strcat('Min number of flagged scans across subjects:', [' ' num2str(min_flagged)],' (',num2str(round(min_flagged_prc,1),'%.1f'),'%)')};
    
    set(FD_MW_info_box,'String',stat_str)
    set(FD_MW_LB1,'String',lb1_str);
end

% Save group statistics & user-specified tmfc_denoise settings
%--------------------------------------------------------------------------
function save_data(~,~)
    if isempty(SPM_paths)
        [filename, pathname] = uiputfile('*.mat', 'Save FD group statistics');
        if isequal(filename,0) || isequal(pathname,0)
            warning('FD group statistics not saved. File name or path not selected.'); 
        else
            fullpath = fullfile(pathname, filename);
            save(fullpath,'FD','FDthr','flagged','N_25prc','N_50prc','N_75prc','N_95prc', ...
                'mean_flagged','sd_flagged','max_flagged','min_flagged');
            fprintf('FD group statistics saved: %s\n', fullpath);
        end
    else   
        [filename, pathname] = uiputfile('*.mat', 'Save FD group statistics and TMFC denoise settings');
        if isequal(filename,0) || isequal(pathname,0)
            warning('FD group statistics not saved. File name or path not selected.'); 
        else
            fullpath = fullfile(pathname, filename);
            denoising_settings.SPM_paths = SPM_paths;
            denoising_settings.subject_paths = subject_paths;
            denoising_settings.options = options;
            denoising_settings.struct_paths = struct_paths;
            denoising_settings.funct_paths = funct_paths; 
            save(fullpath,'denoising_settings','FD','FDthr','flagged', ...
                'N_25prc','N_50prc','N_75prc','N_95prc', ...
                'mean_flagged','sd_flagged','max_flagged','min_flagged');
            fprintf('FD group statistics and TMFC denoise settings saved: %s\n', fullpath);
        end
    end
end

uiwait();
end