# DA-ACh_Dualimaging_CRF
% Close all figures, clear variables
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% Define paths and data to analyze

path2data = '/Users/kellymartyniuk/Desktop/Kellendonk Lab/Projects/Dual imaging/Analysis/Raw data/Eticlopride/'; %Location of the dBASELINE_PERata
path2savefolder = '/Users/kellymartyniuk/Desktop/Kellendonk Lab/Projects/Dual imaging/Analysis/Anazyled data/Eticlopride/'; %Path to save

%% Define the list of mice to analyzed

% Determine it by the folders included in the path2data
mice_list = dir(path2data);
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0
        mice_list(o) = [];
    else
        if strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);

%% Define sessions to exclude 
% For the moment we have to add them manually until we define an automatic 
% rejection method. A cell for each animal with the names of the sessions
Sessions2exclude{1} = []; % [] for adding all
Sessions2exclude{2} = [];
Sessions2exclude{3} = [];

%% Setup the variables for the data you want to extract

IdChannel = {'Grab','dLight'}; % Names of the channels. First = channel A; second channel B
REF_EPOC = 'PC0/'; % event store name. This holds behavioral codes that are
min2remove = 1; % Minutes to remove from the beginning of the recording
TRANGE = [-5 10]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-5 0]; % baseline period within our window
ARTIFACT = Inf; % optionally set an artifact rejection level
THRESHOLD.dLight.peak = 2; % Number of standar deviations to use for the thresholding
THRESHOLD.dLight.dip = 2; % Number of standar deviations to use for the thresholding
THRESHOLD.Grab.dip = 1; % Number of standar deviations to use for the thresholding
THRESHOLD.Grab.peak = 2; % Number of standar deviations to use for the thresholding
Max_lag = 2; % Max time lag in seconds to test for correlation
corr_type = 'Pearson'; % Pearson or Spearman
show_plot = 0; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
limits2plot.dLight = []; %If empty, automatically adjusted to the data. 
limits2plot.Grab = []; %If empty, automatically adjusted to the data.
reanalysis = 1; % If 1, the code runs for sessions already analyze. If 0, session excluded if analysed
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data

%% Loop for all the mice
for m = 1:Nmice
    
%     Define the path to the mouse and find the folders to analyze:
    path2mouse = [path2data,mice_list(m).name,'/'];
    sessions = dir(path2mouse);
    for o = length(sessions):-1:1
        if contains(sessions(o).name,mice_list(m).name) == 0
            sessions(o)=[];
        end
    end
    
    % Create a folder to save the data for this mouse and define the path
    if exist([path2savefolder,mice_list(m).name],'dir') == 0
        mkdir([path2savefolder,mice_list(m).name])
    end
    path2save_mouse = [path2savefolder,mice_list(m).name,'/'];
    
    behav_data2load = dir(path2mouse);
    for o = length(behav_data2load):-1:1
        if behav_data2load(o).isdir == 1
            behav_data2load(o) = [];
        end
    end
    
    %% Loop for all the sessions for the mouse
    for s = 1:length(sessions)
        IdChannel = {'Grab','dLight'};
        % Define the path to the session and create folder to save if needed:
        PATH2SESSION = [path2mouse,sessions(s).name];
        if exist([path2save_mouse,sessions(s).name],'dir') == 0
            mkdir([path2save_mouse,sessions(s).name])
        end
        PATH2SAVE = [path2save_mouse,sessions(s).name,'/'];
        
        % Check if results are already saved for this session 
        done = exist([PATH2SAVE,'ITI_analysis.mat'],'file');
        if done == 0 || reanalysis == 1
            if exist([PATH2SAVE,'figures'],'dir') == 0
                mkdir([PATH2SAVE,'figures'])
            end
            
            %% Read the data
            % Now read the specified data from our block into a Matlab structure.
            data = TDTbin2mat_MAC(PATH2SESSION, 'TYPE', {'epocs', 'scalars', 'streams'});
            
            % Extract the names of the STREAMS in data (A405A or X05A)
            STREAM_STORE1 = cell(1,2);
            STREAM_STORE2 = cell(1,2);
            stream_Names = fieldnames(data.streams);
            for o = 1:length(stream_Names)
                if contains(stream_Names{o},'05A')
                    STREAM_STORE1{1} = stream_Names{o};
                elseif contains(stream_Names{o},'05B')
                    STREAM_STORE1{2} = stream_Names{o};
                elseif contains(stream_Names{o},'65A')
                    STREAM_STORE2{1} = stream_Names{o};
                elseif contains(stream_Names{o},'65B')
                    STREAM_STORE2{2} = stream_Names{o};
                end
            end
            
            % Select the trials based on timestamps and TRANGE (OLD VERSION)
            % data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);
            
            % Get the trig times of the stimulus
            stim_times = [data.epocs.PC0_.onset data.epocs.PC0_.offset];
            
            %% Load the operant box information to align the time:
            % Identify the data from this session (only )
            data2load = [];
            idx = strfind(sessions(s).name,'_');
            session = sessions(s).name(idx(3)+1:end);
            for ii = 1:length(behav_data2load)
                idx = strfind(behav_data2load(ii).name,'_');
                behav_session = behav_data2load(ii).name(idx+1:end-4);
                if strfind(session,behav_session) == 1
                    data2load = ii;
                end
            end
            
            for ii = 1
                if ~isempty(data2load)
                    load([path2mouse,behav_data2load(data2load).name],'Var','rawlist')
                    
                    for nm = 1:length(rawlist.name)
                        if contains(rawlist.name{nm},mice_list(m).name(5:end)) == 1
                            behav_data = nm;
                        end
                    end
                    
                    L_Lever_On_Rows = find(Var{behav_data}(:,2) == 27);
                    R_Lever_On_Rows = find(Var{behav_data}(:,2) == 28);
                    Lever_On_Rows = sort([L_Lever_On_Rows;R_Lever_On_Rows]);
                    t_LeverOn = Var{behav_data}(Lever_On_Rows,1);
                    t_all_stim = t_LeverOn;
                    
                    if size(t_all_stim,1)~= size(stim_times,1)
                        if length(t_all_stim) > length(stim_times)
                            n_errors = length(t_all_stim) - length(stim_times);
                            idx_error = ones(n_errors,1)*nan;
                            
                            tmp1 = stim_times(:,1)-stim_times(1,1);
                            tmp1 = tmp1.\tmp1(end);
                            tmp2 = t_all_stim - t_all_stim(1);
                            tmp2 = tmp2.\tmp2(end);
                            ooo = 0;
                            for oo = 1:length(t_all_stim)
                                temp = tmp1(oo - ooo) - tmp2(oo);
                                if temp >= 0.01
                                    ooo = ooo + 1;
                                    idx_error(ooo) = oo;
                                    if ooo == n_errors
                                        break
                                    end
                                end
                            end
                            t_all_stim(idx_error) = [];
                        else
                            n_errors = length(stim_times) - length(t_all_stim);
                            idx_error = ones(n_errors,1)*nan;
                            
                            tmp1 = stim_times(:,1)-stim_times(1,1);
                            tmp1 = tmp1.\tmp1(end);
                            tmp2 = t_all_stim - t_all_stim(1);
                            tmp2 = tmp2.\tmp2(end);
                            ooo = 0;
                            for oo = 1:length(stim_times)
                                temp = tmp2(oo - ooo) - tmp1(oo);
                                if temp >= 0.01
                                    ooo = ooo + 1;
                                    idx_error(ooo) = oo;
                                    if ooo == n_errors
                                        break
                                    end
                                end
                            end
                            stim_times(idx_error,:) = [];
                        end
                    end
                    
                    % Obtain the linear function for the increasing mismatch
                    % between the fiber photometry and the operant boxes systems
                    x = t_all_stim(:,1);
                    y = stim_times(:,1) - t_all_stim(:,1);
                    P = polyfit(x,y,1);
                    t_diff = P(1)*Var{behav_data}(:,1) + P(2);
                    
                    % Adjust the time events from the operant box to the corrected
                    % time aligned to the fiber photometry
                    Adjusted_Var = Var{behav_data};
                    Adjusted_Var(:,1) = Adjusted_Var(:,1) + t_diff;
                    
                    % Get the new timing for the head entries and stims
                    L_Lever_press_Rows = find(Var{behav_data}(:,2) == 1015);
                    R_Lever_press_Rows = find(Var{behav_data}(:,2) == 1016);
                    Lever_press_Rows = sort([L_Lever_press_Rows;R_Lever_press_Rows]);
                    L_Lever_Off_Rows = find(Var{behav_data}(:,2) == 29);
                    R_Lever_Off_Rows = find(Var{behav_data}(:,2) == 30);
                    Lever_Off_Rows = sort([L_Lever_Off_Rows;R_Lever_Off_Rows]);
                    Dipper_On_Rows = find(Var{behav_data}(:,2) == 25);
                    Dipper_Off_Rows = find(Var{behav_data}(:,2) == 26);
                    Head_entries_Rows = find(Adjusted_Var(:,2) == 1011);
                    if length(Lever_On_Rows) == length(Lever_Off_Rows)
                        t_Lever = [Adjusted_Var(Lever_On_Rows,1) Adjusted_Var(Lever_Off_Rows,1)];
                    elseif length(Lever_On_Rows) > length(Lever_Off_Rows)
                        t_Lever = ones(length(Lever_On_Rows),2)*nan;
                        t_Lever(:,1) = Adjusted_Var(Lever_On_Rows,1);
                        t_Lever(1:length(Lever_Off_Rows),2) = Adjusted_Var(Lever_Off_Rows,1);
                    end
                    t_LeverPress = Adjusted_Var(Lever_press_Rows,1);
                    if length(Dipper_On_Rows) == length(Dipper_Off_Rows)
                        t_Dipper = [Adjusted_Var(Dipper_On_Rows,1) Adjusted_Var(Dipper_Off_Rows,1)];
                    elseif length(Dipper_On_Rows) > length(Dipper_Off_Rows)
                        t_Dipper = ones(length(Dipper_On_Rows),2)*nan;
                        t_Dipper(:,1) = Adjusted_Var(Dipper_On_Rows,1);
                        t_Dipper(1:length(Dipper_Off_Rows),2) = Adjusted_Var(Dipper_Off_Rows,1);
                    end
                    [~,idx] = setdiff(t_Dipper(:,2),t_Lever(:,2));
                    if ~isempty(idx)
                        t_Dipper(idx,:) = [];
                    end
                    t_head_entry = Adjusted_Var(Head_entries_Rows,1);
                    
                    % Reject head entries happenning within half a second
                    A = diff(t_head_entry);
                    B = find(A < 0.5);
                    t_head_entry(B+1) = [];
                else
                    t_Lever = [data.epocs.PC0_.onset data.epocs.PC0_.offset];
                    t_LeverPress = [];
                    t_Dipper = [data.epocs.PC0_.offset-5 data.epocs.PC0_.offset];
                    t_head_entry = [];
                end
            end
            
            %Adjust the lengths of both channels
            minLength = min([length(data.streams.(STREAM_STORE1{1}).data) ...
                length(data.streams.(STREAM_STORE1{2}).data) ...
                length(data.streams.(STREAM_STORE2{1}).data)...
                length(data.streams.(STREAM_STORE2{2}).data)]);
            
            for channel = 1:length(STREAM_STORE1)
                if length(data.streams.(STREAM_STORE1{channel}).data) > minLength
                    data.streams.(STREAM_STORE1{channel}).data...
                        (minLength+1:length(data.streams.(STREAM_STORE1{channel}).data)) = [];
                    %(1:length(data.streams.(STREAM_STORE1{channel}).data)-minLength) = [];
                    
                end
                if length(data.streams.(STREAM_STORE2{channel}).data) > minLength
                    data.streams.(STREAM_STORE2{channel}).data...
                        (minLength+1:length(data.streams.(STREAM_STORE2{channel}).data)) = [];
                    %                     (1:length(data.streams.(STREAM_STORE2{channel}).data)-minLength) = [];
                end
            end
            
            % Gets time vector for our data
            Fs = data.streams.(STREAM_STORE1{1}).fs;
            Max_idx = minLength;
            dt = 1/Fs;
            time = 0:dt:dt*(Max_idx-1);
            
            dummie = 1:length(time);
            dummie = dummie(time > min2remove*60);
            idx_remove = dummie(1);
            clear dummie
            
            time_dwn = downsample(time(idx_remove:end),10);
            N = length(time_dwn);
            
            clear time
            
            %% Get the dFF for both sensors
            % Filter out the slow fluctuation
            ftype = 'high';
            n = 2;
            Wn = 0.05/((Fs/10)/2);
            [a,b] = butter(n,Wn,ftype);
            
            for channel = 1:length(STREAM_STORE1)
                
                % downsample 10x
                data405 = data.streams.(STREAM_STORE1{channel}).data(idx_remove:end);
                data405 = downsample(data405,10);
                
                data465 = data.streams.(STREAM_STORE2{channel}).data(idx_remove:end);
                data465 = downsample(data465,10);
                
                %Normalization of the 465 using the 405
                bls = polyfit(data405, data465, 1);
                Y_fit = bls(1) .* data405 + bls(2);
                Y_dF = data465 - Y_fit;
                
                %dFF using 405 fit as baseline
                dFF.(IdChannel{channel}).raw = 100*(Y_dF)./Y_fit;
                dFF.(IdChannel{channel}).zscore = zscore(dFF.(IdChannel{channel}).raw);
                dFF.(IdChannel{channel}).filt = filtfilt(a,b,double(dFF.(IdChannel{channel}).raw));
                dFF.(IdChannel{channel}).zscorefilt = filtfilt(a,b,double(dFF.(IdChannel{channel}).zscore));
                
            end
            
            %% Visualize the whole session:
            plotAllSession = 1; % If 0 don't plot
            for ooo = 1
                if plotAllSession
                    figure
                    subplot(2,1,1)
                    plot(time_dwn,dFF.dLight.filt,'b')
                    hold on
                    plot(time_dwn,dFF.Grab.filt,'r')
                    for o = 1:size(t_Lever,1)
                        if o == 1
                            if ~isempty(t_Lever)
                                plot([t_Lever(o,1) t_Lever(o,1)],[-4 12],'m')
                                plot([t_Lever(o,2) t_Lever(o,2)],[-4 12],'k')
                            end
                            if ~isempty(t_Dipper)
                                plot([t_Dipper(o,1) t_Dipper(o,1)],[-4 12],'c')
                                plot([t_Dipper(o,2) t_Dipper(o,2)],[-4 12],'b')
                            end
                        else
                            if ~isempty(t_Lever)
                                plot([t_Lever(o,1) t_Lever(o,1)],[-4 12],'m','HandleVisibility','off')
                                plot([t_Lever(o,2) t_Lever(o,2)],[-4 12],'k','HandleVisibility','off')
                            end
                            if ~isempty(t_Dipper)
                                if size(t_Dipper,1) >= o
                                    plot([t_Dipper(o,1) t_Dipper(o,1)],[-4 12],'c','HandleVisibility','off')
                                    plot([t_Dipper(o,2) t_Dipper(o,2)],[-4 12],'b','HandleVisibility','off')
                                end
                            end
                        end
                    end
                    if ~isempty(data2load)
                         if ~isempty(t_LeverPress)
                        plot(t_LeverPress,ones(length(t_LeverPress),1)*10,'*g')
                         end 
                        plot(t_head_entry,ones(length(t_head_entry),1)*10,'*k')
                    end
                    xlim([time_dwn(1) time_dwn(end)])
                    xlabel('Time (s)')
                    ylim([-10 30])
                    ylabel('Signal size (A.U.)')
                    title('Raw')
                    if ~isempty(data2load)
                        legend('dLight','Grab','Lever On','Lever Off','Dipper On','Dipper Off','Lever Press','Head entries','Location','northoutside','NumColumns',6)
                    else
                        legend('dLight','Grab','Lever On','Lever Off','Dipper On','Dipper Off','Location','northoutside','NumColumns',4)
                    end
                    subplot(2,1,2)
                    plot(time_dwn,dFF.dLight.zscorefilt,'b')
                    hold on
                    plot(time_dwn,dFF.Grab.zscorefilt,'r')
                    for o = 1:size(t_Lever,1)
                        if o == 1
                            if ~isempty(t_Lever)
                                plot([t_Lever(o,1) t_Lever(o,1)],[-4 12],'m')
                                plot([t_Lever(o,2) t_Lever(o,2)],[-4 12],'k')
                            end
                            if ~isempty(t_Dipper)
                                plot([t_Dipper(o,1) t_Dipper(o,1)],[-4 12],'c')
                                plot([t_Dipper(o,2) t_Dipper(o,2)],[-4 12],'b')
                            end
                        else
                            if ~isempty(t_Lever)
                                plot([t_Lever(o,1) t_Lever(o,1)],[-4 12],'m','HandleVisibility','off')
                                plot([t_Lever(o,2) t_Lever(o,2)],[-4 12],'k','HandleVisibility','off')
                            end
                            if size(t_Dipper,1) >= o
                                if ~isempty(t_Dipper)
                                    plot([t_Dipper(o,1) t_Dipper(o,1)],[-4 12],'c','HandleVisibility','off')
                                    plot([t_Dipper(o,2) t_Dipper(o,2)],[-4 12],'b','HandleVisibility','off')
                                end
                            end
                        end
                    end
                    if ~isempty(data2load)
                        if ~isempty(t_LeverPress)
                            plot(t_LeverPress,ones(length(t_LeverPress),1)*10,'*g')
                        end
                        plot(t_head_entry,ones(length(t_head_entry),1)*10,'*k')
                    end
                    xlim([time_dwn(1) time_dwn(end)])
                    xlabel('Time (s)')
                    ylim([-5 10])
                    ylabel('Signal size (A.U.)')
                    title('Z-scored')
                    if ~isempty(data2load)
                        legend('Zscore dLight','Zscore Grab','Lever On','Lever Off','Dipper On','Dipper Off','Lever Press','Head entries','Location','northoutside','NumColumns',6)
                    else
                        legend('Zscore dLight','Zscore Grab','Lever On','Lever Off','Dipper On','Dipper Off','Location','northoutside','NumColumns',4)
                    end
                    if save_plot == 1
                        saveas(gcf,[PATH2SAVE,'figures/Full session.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures/Full session.fig'])
                    end
                end
            end
            
            %% Get moving window correlation
            winsize = 2;
            dt = 0.5;
            
            temp_t = time_dwn - time_dwn(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= winsize);
            idx_winsize = dummie(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= dt);
            idx_dt = dummie(1);
            clear dummie temp_t
            
                        
            % Centered in the time window
            half_win = round(idx_winsize\2);
            half_dt = round(idx_dt\2);
            tmpix = half_win - half_dt;
            
            init_ix = 1:idx_dt:(length(time_dwn)-idx_winsize);
            end_ix = idx_dt:idx_dt:(length(time_dwn)-idx_winsize+idx_dt);
            end_ix(end) = length(time_dwn);
            
            signal1 = dFF.(IdChannel{1}).filt;
            signal2 = dFF.(IdChannel{2}).filt;
            moving_corr = ones(1,length(time_dwn))*nan;
            for o = 1:length(init_ix)
                x = signal1(init_ix(o):init_ix(o)+idx_winsize);
                y = signal2(init_ix(o):init_ix(o)+idx_winsize);
                [r,~] = corr(x',y','Type',corr_type);
                if o == 1
                    moving_corr(1:init_ix(o)+idx_winsize) = r;
                elseif o == length(init_ix)
                    moving_corr((init_ix(o)+tmpix):end) = r;
                else
                    moving_corr((init_ix(o)+tmpix):init_ix(o)+idx_winsize) = r;
                end
            end
            clear idx_winsize idx_dt half_win half_dt tmpix init_ix end_ix signal1 signal2
            
            %% Get overall correlation 
            %%% Using raw data (filt version)
            signal1 = dFF.dLight.filt;
            signal2 = dFF.Grab.filt;
            save_name = 'Raw overall correlation with Grab lag';
            
            % Compute the correlation when lagging Grab agaist dLight
            [lag2plot,Ovrl_corr.raw] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            clear signal1 signal2
            
            %%% Using zscore data (filt version)
            signal1 = dFF.dLight.zscorefilt;
            signal2 = dFF.Grab.zscorefilt;
            save_name = 'Zscore overall correlation with Grab lag';
            
            % Compute the correlation when lagging Grab agaist dLight
            [~,Ovrl_corr.zscore] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            
            clear signal1 signal2
            close all
            
            %% Trials analysis
            % Define time index for the trial period
            temp_t = time_dwn - time_dwn(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= abs(TRANGE(1)));
            idx_Init = dummie(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= abs(TRANGE(2)));
            idx_End = dummie(1);
            
            % Define time index for baseline correction
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= abs(BASELINE_PER(1)));
            baseline_idx1 = dummie(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= abs(BASELINE_PER(2)));
            baseline_idx2 = dummie(1);
            
            n = idx_Init + idx_End; % Length of each trial
            t_trials = temp_t(1:n) - temp_t(idx_Init);
            
            % Find first head entry after the dipper on to align for the
            % reward
            latency.leverExtention2press = ones(size(t_Lever,1),1)*nan;
            latency.leverExtention2pressMore2s = ones(size(t_Lever,1),1)*nan;
            latency.leverExtention2pressMore1s = ones(size(t_Lever,1),1)*nan;
            latency.leverExtention2reward = ones(size(t_Lever,1),1)*nan;
            latency.leverPress2reward = ones(size(t_Lever,1),1)*nan;
            t_reward = ones(size(t_Dipper,1),1)*nan;
            

            
            for o = 1:size(t_Lever,1)
                if o <= size(t_Dipper,1)
                    latency.leverExtention2press(o) = t_Dipper(o,1) - t_Lever(o,1);
                    if ~isempty(data2load)
                        tmp_idx = t_head_entry > t_Dipper(o,1) & t_head_entry < t_Dipper(o,2);
                        tmp_time = t_head_entry(tmp_idx);
                        if ~isempty(tmp_time)
                            t_reward(o) = tmp_time(1);
                            latency.leverExtention2reward(o) = tmp_time(1) - t_Lever(o,1);
                            latency.leverPress2reward(o) = tmp_time(1) - t_Dipper(o,1);
                        end
                    end
                end
            end
            t_reward(isnan(t_reward)) = [];
            latency.leverExtention2press(isnan(latency.leverExtention2press))= [];
            latency.leverExtention2pressMore2s= find(latency.leverExtention2press > 2.0);
            latency.leverExtention2pressMore1s= find(latency.leverExtention2press > 1.0);
            [~,latency_sorted]= sort (latency.leverExtention2press);
            mean_LPlatency= mean(latency.leverExtention2press);
            
            %% Compute the aligned results of the trials
            Stim_data.LeverExtension.idx = ones(size(t_Lever,1),n)*nan;
            Stim_data.LeverExtension.corr = ones(size(t_Lever,1),n)*nan;
            Stim_data.LeverExtension.head_entries.idx = cell(size(t_Lever,1),1);
            Stim_data.LeverExtension.head_entries.times = cell(size(t_Lever,1),1);
            Stim_data.LeverExtension.LeverPress.idx = cell(size(t_Lever,1),1);
            Stim_data.LeverExtension.LeverPress.times = cell(size(t_Lever,1),1);
            
            Stim_data.DipperOn.idx = ones(size(t_Dipper,1),n)*nan;
            Stim_data.DipperOn.corr = ones(size(t_Dipper,1),n)*nan;
            Stim_data.DipperOn.head_entries.idx = cell(size(t_Dipper,1),1);
            Stim_data.DipperOn.head_entries.times = cell(size(t_Dipper,1),1);
            Stim_data.DipperOn.LeverPress.idx = cell(size(t_Dipper,1),1);
            Stim_data.DipperOn.LeverPress.idx = cell(size(t_Dipper,1),1);
            Stim_data.DipperOn.LeverPress.times = cell(size(t_Dipper,1),1);
            if ~isempty(data2load)
                Stim_data.Reward.idx = ones(size(t_reward,1),n)*nan;
                Stim_data.Reward.corr = ones(size(t_reward,1),n)*nan;
                Stim_data.Reward.head_entries.idx = cell(size(t_reward,1),1);
                Stim_data.Reward.head_entries.times = cell(size(t_reward,1),1);
                Stim_data.Reward.LeverPress.idx = cell(size(t_reward,1),1);
                Stim_data.Reward.LeverPress.times = cell(size(t_reward,1),1);
            end
            for i = 1:length(IdChannel)
                Stim_data.LeverExtension.dFF.(IdChannel{i}).raw = ones(size(t_Lever,1),n)*nan;
                Stim_data.LeverExtension.dFF.(IdChannel{i}).baseline_corrected = ones(size(t_Lever,1),n)*nan;
                Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore = ones(size(t_Lever,1),n)*nan;
                Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore_baseline_corrected = ones(size(t_Lever,1),n)*nan;
                
                Stim_data.DipperOn.dFF.(IdChannel{i}).raw = ones(size(t_Dipper,1),n)*nan;
                Stim_data.DipperOn.dFF.(IdChannel{i}).baseline_corrected = ones(size(t_Dipper,1),n)*nan;
                Stim_data.DipperOn.dFF.(IdChannel{i}).zscore = ones(size(t_Dipper,1),n)*nan;
                Stim_data.DipperOn.dFF.(IdChannel{i}).zscore_baseline_corrected = ones(size(t_Dipper,1),n)*nan;
                
                if ~isempty(data2load)
                    Stim_data.Reward.dFF.(IdChannel{i}).raw = ones(size(t_reward,1),n)*nan;
                    Stim_data.Reward.dFF.(IdChannel{i}).baseline_corrected = ones(size(t_reward,1),n)*nan;
                    Stim_data.Reward.dFF.(IdChannel{i}).zscore = ones(size(t_reward,1),n)*nan;
                    Stim_data.Reward.dFF.(IdChannel{i}).zscore_baseline_corrected = ones(size(t_reward,1),n)*nan;
                end
            end
            
            % Lever Extension
            for o = 1:size(t_Lever,1)
                ix = find(abs(time_dwn-t_Lever(o,1)) == min(abs(time_dwn-t_Lever(o,1))));
                if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                    tmp = ix - (idx_Init-1):ix + idx_End;
                    Stim_data.LeverExtension.idx(o,:) = tmp;
                    for i = 1:length(IdChannel)
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).raw(o,:) = ...
                            dFF.(IdChannel{i}).filt(tmp);
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore(o,:) = ...
                            dFF.(IdChannel{i}).zscorefilt(tmp);
                    end
                    Stim_data.LeverExtension.corr(o,:) = moving_corr(tmp);
                elseif ix+idx_End > length(time_dwn)
                    tmp = ix - (idx_Init-1):length(time_dwn);
                    Stim_data.LeverExtension.idx(o,1:length(tmp)) = tmp;
                    for i = 1:length(IdChannel)
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).raw(o,1:length(tmp)) = ...
                            dFF.(IdChannel{i}).filt(tmp);
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore(o,1:length(tmp)) = ...
                            dFF.(IdChannel{i}).zscorefilt(tmp);
                    end
                    Stim_data.LeverExtension.corr(o,1:length(tmp)) = moving_corr(tmp);
                elseif ix-(idx_Init-1) < 1
                    tmp = 1:(ix + idx_End);
                    Stim_data.LeverExtension.idx(o,(n-length(tmp)+1):end) = tmp;
                    for i = 1:length(IdChannel)
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).raw(o,(n-length(tmp)+1):end) = ...
                            dFF.(IdChannel{i}).filt(tmp);
                        Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore(o,(n-length(tmp)+1):end) = ...
                            dFF.(IdChannel{i}).zscorefilt(tmp);
                    end
                    Stim_data.LeverExtension.corr(o,(n-length(tmp)+1):end) = moving_corr(tmp);
                end
                if ~isempty(data2load)
                    dummie = 1:length(t_head_entry);
                    Stim_data.LeverExtension.head_entries.idx{o} = dummie...
                        (t_head_entry >= time_dwn(tmp(1)) & t_head_entry <= time_dwn(tmp(end)));
                    Stim_data.LeverExtension.head_entries.times{o} = ...
                        t_head_entry(Stim_data.LeverExtension.head_entries.idx{o})-t_Lever(o,1);
                    
                    dummie = 1:length(t_LeverPress);
                    Stim_data.LeverExtension.LeverPress.idx{o} = dummie...
                        (t_LeverPress >= time_dwn(tmp(1)) & t_LeverPress <= time_dwn(tmp(end)));
                    Stim_data.LeverExtension.LeverPress.times{o} = ...
                        t_LeverPress(Stim_data.LeverExtension.LeverPress.idx{o})-t_Lever(o);
                end
            end
            for i = 1:length(IdChannel)
                Stim_data.LeverExtension.dFF.(IdChannel{i}).baseline_corrected = ...
                    Stim_data.LeverExtension.dFF.(IdChannel{i}).raw - nanmedian...
                    (Stim_data.LeverExtension.dFF.(IdChannel{i}).raw...
                    (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore_baseline_corrected = ...
                    Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore - nanmedian...
                    (Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore...
                    (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
            end
            
            %% Find trials with lever press latency >2s (alinged to lever extension)
            StimData_2sLat_GACh= Stim_data.LeverExtension.dFF.Grab.baseline_corrected([latency.leverExtention2pressMore2s], :);
            StimData_2sLat_GACh_mean= mean(StimData_2sLat_GACh,'omitnan');
            StimData_2sLat_DA= Stim_data.LeverExtension.dFF.dLight.baseline_corrected ([latency.leverExtention2pressMore2s], :);
            StimData_2sLat_DA_mean = mean(StimData_2sLat_DA,'omitnan');
            
             % All trials regardless of press latency (alinged to lever
            % extension)
            StimData_GACh= Stim_data.LeverExtension.dFF.Grab.baseline_corrected;
            StimData_GACh_mean= mean(StimData_GACh,'omitnan');
            StimData_DA= Stim_data.LeverExtension.dFF.dLight.baseline_corrected ;
            StimData_DA_mean= mean(StimData_DA,'omitnan');
           
            
            % Create heatmap for trials sorted by lever press latency > 2 sec (alinged to lever extension)
            StimData_sorted_LPL_GACh= Stim_data.LeverExtension.dFF.Grab.baseline_corrected ([latency_sorted],:);
            StimData_sorted_LPL_DA= Stim_data.LeverExtension.dFF.dLight.baseline_corrected([latency_sorted],:);
            
            %Creat heatmap for GACh
            figure;
            imagesc(t_trials, 1, StimData_sorted_LPL_GACh);
            colormap('jet'); % c1 = colorbar;
            colorbar;
            line([0 0], [0 60], 'Color', 'k', 'LineStyle','-', 'LineWidth', 2)
            xlim([-2 10]);
            save_name = 'Heat Map Sorted by LPL';
            saveas(gcf,[PATH2SAVE,'/',save_name,'Heat Map Sorted by LPL__ACh.fig'])
            saveas(gcf,[PATH2SAVE,'/',save_name,'Heat Map Sorted by LPL__ACh.jpg'])
            
            %Creat heatmap for DA
            figure;
            imagesc(t_trials, 1, StimData_sorted_LPL_DA);
            colormap('jet'); % c1 = colorbar;
            colorbar;
            line([0 0], [0 60], 'Color', 'k', 'LineStyle','-', 'LineWidth', 2)
            xlim([-2 10]);
            save_name = 'Heat Map Sorted by LPL';
            saveas(gcf,[PATH2SAVE,'/',save_name,'Heat Map Sorted by LPL__DA.fig'])
            saveas(gcf,[PATH2SAVE,'/',save_name,'Heat Map Sorted by LPL__DA.jpg'])

            
            
            %%
            % Dipper On (FIRST LEVER PRESS)
            if ~isempty(t_Dipper)
                for o = 1:size(t_Dipper,1)
                    ix = find(abs(time_dwn-t_Dipper(o,1)) == min(abs(time_dwn-t_Dipper(o,1))));
                    if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                        tmp = ix - (idx_Init-1):ix + idx_End;
                        Stim_data.DipperOn.idx(o,:) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.DipperOn.dFF.(IdChannel{i}).raw(o,:) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.DipperOn.dFF.(IdChannel{i}).zscore(o,:) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.DipperOn.corr(o,:) = moving_corr(tmp);
                    elseif ix+idx_End > length(time_dwn)
                        tmp = ix - (idx_Init-1):length(time_dwn);
                        Stim_data.DipperOn.idx(o,1:length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.DipperOn.dFF.(IdChannel{i}).raw(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.DipperOn.dFF.(IdChannel{i}).zscore(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.DipperOn.corr(o,1:length(tmp)) = moving_corr(tmp);
                    elseif ix-(idx_Init-1) < 1
                        tmp = 1:ix + idx_End;
                        Stim_data.DipperOn.idx(o,(n-length(tmp)+1):end) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.DipperOn.dFF.(IdChannel{i}).raw(o,(n-length(tmp)+1):end) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.DipperOn.dFF.(IdChannel{i}).zscore(o,(n-length(tmp)+1):end) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.DipperOn.corr(o,(n-length(tmp)+1):end) = moving_corr(tmp);
                    end
                    if ~isempty(data2load)
                        dummie = 1:length(t_head_entry);
                        Stim_data.DipperOn.head_entries.idx{o} = dummie...
                            (t_head_entry >= time_dwn(tmp(1)) & t_head_entry <= time_dwn(tmp(end)));
                        Stim_data.DipperOn.head_entries.times{o} = ...
                            t_head_entry(Stim_data.DipperOn.head_entries.idx{o})-t_Dipper(o,1);
                        
                        dummie = 1:length(t_LeverPress);
                        Stim_data.DipperOn.LeverPress.idx{o} = dummie...
                            (t_LeverPress >= time_dwn(tmp(1)) & t_LeverPress <= time_dwn(tmp(end)));
                        Stim_data.DipperOn.LeverPress.times{o} = ...
                            t_LeverPress(Stim_data.DipperOn.LeverPress.idx{o})-t_Dipper(o);
                    end
                end
                for i = 1:length(IdChannel)         
                    Stim_data.DipperOn.dFF.(IdChannel{i}).baseline_corrected = ...
                        Stim_data.DipperOn.dFF.(IdChannel{i}).raw - nanmedian...
                        (Stim_data.DipperOn.dFF.(IdChannel{i}).raw...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                    Stim_data.DipperOn.dFF.(IdChannel{i}).zscore_baseline_corrected = ...
                        Stim_data.DipperOn.dFF.(IdChannel{i}).zscore - nanmedian...
                        (Stim_data.DipperOn.dFF.(IdChannel{i}).zscore...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                end
            end
            %%  % Find trials with lever press latency >2s (alinged to first lever press)
            StimData_LP_2sLat_GACh= Stim_data.DipperOn.dFF.Grab.baseline_corrected([latency.leverExtention2pressMore2s], :);
            StimData_LP_2sLat_GACh_mean= mean(StimData_LP_2sLat_GACh,'omitnan'); 
            StimData_LP_2sLat_DA= Stim_data.DipperOn.dFF.dLight.baseline_corrected ([latency.leverExtention2pressMore2s], :);
            StimData_LP_2sLat_DA_mean = mean(StimData_LP_2sLat_DA,'omitnan');
            
            % All trials regardless of press latency (alinged to lever press)
            StimData_LP_GACh= Stim_data.DipperOn.dFF.Grab.baseline_corrected;
            StimData_LP_GACh_mean= mean(StimData_LP_GACh,'omitnan');
            StimData_LP_DA= Stim_data.DipperOn.dFF.dLight.baseline_corrected ;
            StimData_LP_DA_mean= mean(StimData_LP_DA,'omitnan');
                
            %% 
            % Reward
            if ~isempty(data2load) && ~isempty(t_reward)
                for o = 1:size(t_reward,1)
                    ix = find(abs(time_dwn-t_reward(o,1)) == min(abs(time_dwn-t_reward(o,1))));
                    if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                        tmp = ix - (idx_Init-1):ix + idx_End;
                        Stim_data.Reward.idx(o,:) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.Reward.dFF.(IdChannel{i}).raw(o,:) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.Reward.dFF.(IdChannel{i}).zscore(o,:) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.Reward.corr(o,:) = moving_corr(tmp);
                    elseif ix+idx_End > length(time_dwn)
                        tmp = ix - (idx_Init-1):length(time_dwn);
                        Stim_data.Reward.idx(o,1:length(tmp)) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.Reward.dFF.(IdChannel{i}).raw(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.Reward.dFF.(IdChannel{i}).zscore(o,1:length(tmp)) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.Reward.corr(o,1:length(tmp)) = moving_corr(tmp);
                    elseif ix-(idx_Init-1) < 1
                        tmp = 1:ix + idx_End;
                        Stim_data.Reward.idx(o,(n-length(tmp)+1):end) = tmp;
                        for i = 1:length(IdChannel)
                            Stim_data.Reward.dFF.(IdChannel{i}).raw(o,(n-length(tmp)+1):end) = ...
                                dFF.(IdChannel{i}).filt(tmp);
                            Stim_data.Reward.dFF.(IdChannel{i}).zscore(o,(n-length(tmp)+1):end) = ...
                                dFF.(IdChannel{i}).zscorefilt(tmp);
                        end
                        Stim_data.Reward.corr(o,(n-length(tmp)+1):end) = moving_corr(tmp);
                    end
                    if ~isempty(data2load)
                        dummie = 1:length(t_head_entry);
                        Stim_data.Reward.head_entries.idx{o} = dummie...
                            (t_head_entry >= time_dwn(tmp(1)) & t_head_entry <= time_dwn(tmp(end)));
                        Stim_data.Reward.head_entries.times{o} = ...
                            t_head_entry(Stim_data.Reward.head_entries.idx{o})-t_reward(o,1);
                        
                        dummie = 1:length(t_LeverPress);
                        Stim_data.Reward.LeverPress.idx{o} = dummie...
                            (t_LeverPress >= time_dwn(tmp(1)) & t_LeverPress <= time_dwn(tmp(end)));
                        Stim_data.Reward.LeverPress.times{o} = ...
                            t_LeverPress(Stim_data.Reward.LeverPress.idx{o})-t_reward(o,1);
                    end
                end
                for i = 1:length(IdChannel)
                    Stim_data.Reward.dFF.(IdChannel{i}).baseline_corrected = ...
                        Stim_data.Reward.dFF.(IdChannel{i}).raw - nanmedian...
                        (Stim_data.Reward.dFF.(IdChannel{i}).raw...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                    Stim_data.Reward.dFF.(IdChannel{i}).zscore_baseline_corrected = ...
                        Stim_data.Reward.dFF.(IdChannel{i}).zscore - nanmedian...
                        (Stim_data.Reward.dFF.(IdChannel{i}).zscore...
                        (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                end
            end
            clear dummie temp_t
            
            % Plot the results
            norm = {'raw','zscore'};
            Condition = fieldnames(Stim_data);%{'LeverExtension','DipperOn','Reward'};
            for cond = 1:length(Condition)
                if sum(~isnan(Stim_data.(Condition{cond}).dFF.(IdChannel{1}).baseline_corrected(:,760))) ~= 0
                    
                    %%% Raw
                    Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).baseline_corrected;
                    Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).baseline_corrected;
                    plot_traces_and_behavior_AT(Signal1,Signal2,t_trials,Stim_data.(Condition{cond}).head_entries.times,...
                        Stim_data.(Condition{cond}).LeverPress.times,Condition{cond},norm{1},IdChannel,TRANGE,limits2plot.(IdChannel{1}),...
                        limits2plot.(IdChannel{2}),PATH2SAVE,show_plot,save_plot,data2load)
                    
                    %%% Z-score
                    Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).zscore_baseline_corrected;
                    Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).zscore_baseline_corrected;
                    plot_traces_and_behavior_AT(Signal1,Signal2,t_trials,Stim_data.(Condition{cond}).head_entries.times,...
                        Stim_data.(Condition{cond}).LeverPress.times,Condition{cond},norm{2},IdChannel,TRANGE,limits2plot.(IdChannel{1}),...
                        limits2plot.(IdChannel{2}),PATH2SAVE,show_plot,save_plot,data2load)
                    
                    %% Plot the average responses and the moving correlation
                    %%% Raw
                    Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).baseline_corrected;
                    Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).baseline_corrected;
                    plot_average_traces_and_correlation_AT(Signal1,Signal2,Stim_data.(Condition{cond}).corr,...
                        t_trials,Condition{cond},norm{1},PATH2SAVE,show_plot,save_plot)
                    %%% Z-score
                    Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).zscore_baseline_corrected;
                    Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).zscore_baseline_corrected;
                    plot_average_traces_and_correlation_AT(Signal1,Signal2,Stim_data.(Condition{cond}).corr,...
                        t_trials,Condition{cond},norm{2},PATH2SAVE,show_plot,save_plot)
                end
            end
% %             
            %% Compute the AUC and the Peaks (lever extension)
            temp_t = time_dwn - time_dwn(1);
              idx_AUC1 = find(t_trials >= 0);
              idx_AUC1= idx_AUC1 (1);
              idx_AUC2= find(t_trials >= 5);
              idx_AUC2= idx_AUC2 (1);
              idx_AUC= idx_AUC1:idx_AUC2;
            clear dummie idx_AUC1 idx_AUC2
            
            for i = 1:length(IdChannel)
                Measurements.LeverExtension.(IdChannel{i}).raw.AUC = cell(size(t_Lever,1),1);
                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks = cell(size(t_Lever,1),1);
                selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.val = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.time = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.val = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.time = ones(size(t_Lever,1),3)*nan;
                
                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC = cell(size(t_Lever,1),1);
                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks = cell(size(t_Lever,1),1);
                selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.val = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.time = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.val = ones(size(t_Lever,1),3)*nan;
                selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.time = ones(size(t_Lever,1),3)*nan;
            end
            % Divide the trace in each sector above or bellow 0. Get the 
            % AUC and peak for each section
            for o = 1:size(t_Lever,1)
                for i = 1:length(IdChannel)
                    %%% Raw
                    tmp = Stim_data.LeverExtension.dFF.(IdChannel{i}).baseline_corrected(o,idx_AUC);
                    A = tmp >= 0;
                    B = find(A == 1);
                    if ~isempty(B)
                        C = diff(B);
                        D = find(C ~= 1);
                        zero_crossing = unique(sort([1 B(1) B(D) B(D+1) B(end) length(tmp)]));
                        tmp_AUC = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.val = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.idx = ones(length(zero_crossing)-1,1)*nan;
                        for oo = 1:length(zero_crossing)-1
                            tmp2 = tmp(zero_crossing(oo)+1:zero_crossing(oo+1));
                            tmp_AUC(oo) = trapz(tmp2);
                            [~,idx] = max(abs(tmp2));
                            tmp_peaks.val(oo) = tmp2(idx);
                            tmp_peaks.idx(oo) = zero_crossing(oo)+idx-1;
                        end
                        tmp2 = temp_t(zero_crossing);
                        Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o} = ...
                            [tmp_AUC diff(tmp2)' tmp2(1:end-1)'];
                        Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o} = ...
                            [tmp_peaks.val temp_t(tmp_peaks.idx)'];
                    else
                        [~,idx] = max(abs(tmp));
                        tmp_peaks.val(oo) = tmp(idx);
                        tmp_peaks.idx(oo) = idx;
                        Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o} = ...
                            [tmp_peaks.val(oo) temp_t(tmp_peaks.idx(oo))'];
                        Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o} = ...
                            [trapz(tmp) temp_t(length(idx_AUC)) 0];
                    end
                    
                    %%% Zscore
                    tmp = Stim_data.LeverExtension.dFF.(IdChannel{i}).zscore_baseline_corrected(o,idx_AUC);
                    A = tmp >= 0;
                    B = find(A == 1);
                    if ~isempty(B)
                        C = diff(B);
                        D = find(C ~= 1);
                        zero_crossing = unique(sort([1 B(1) B(D) B(D+1) B(end) length(tmp)]));
                        tmp_AUC = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.val = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.idx = ones(length(zero_crossing)-1,1)*nan;
                        for oo = 1:length(zero_crossing)-1
                            tmp2 = tmp(zero_crossing(oo)+1:zero_crossing(oo+1));
                            tmp_AUC(oo) = trapz(tmp2);
                            [~,idx] = max(abs(tmp2));
                            tmp_peaks.val(oo) = tmp2(idx);
                            tmp_peaks.idx(oo) = zero_crossing(oo)+idx-1;
                        end
                        tmp2 = temp_t(zero_crossing);
                        Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o} = ...
                            [tmp_AUC diff(tmp2)' tmp2(1:end-1)'];
                        Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o} = ...
                            [tmp_peaks.val temp_t(tmp_peaks.idx)'];
                    else
                        [~,idx] = max(abs(tmp));
                        tmp_peaks.val(oo) = tmp(idx);
                        tmp_peaks.idx(oo) = idx;
                        Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o} = ...
                            [tmp_peaks.val(oo) temp_t(tmp_peaks.idx(oo))'];
                        Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o} = ...
                            [trapz(tmp) temp_t(length(idx_AUC)) 0];
                    end
                end
            end
            % Select the relevant sections: negative-positive-negative for 
            % DLight,positive-negative-positive for GaCh 
            for o = 1:size(t_Lever,1)
                for i = 1:length(IdChannel)
                    %%% Raw
                    if contains(IdChannel{i},'Grab')
                        A = find(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(:,1) < 0 ...
                            & Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(:,3) < 1);
                    else
                        A = find(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(:,1) > 0 ...
                            & Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(:,3) < 1);
                    end
                    if ~isempty(A)
                        [~,idx] = max(abs(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A,1)));
                        A = A(idx);
                        if A > 1 && A < size(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.val(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A-1:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.time(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A-1:A+1,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.val(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A-1:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.time(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A-1:A+1,2);
                        elseif A > 1 && A >= size(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.val(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A-1:A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.time(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A-1:A,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.val(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A-1:A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.time(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A-1:A,2);
                        elseif A == 1 && A < size(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.val(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.time(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A:A+1,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.val(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.time(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A:A+1,2);
                        elseif A == 1 && A >= size(Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.val(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.AUC.time(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.AUC{o}(A,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.val(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).raw.Peaks.time(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).raw.Peaks{o}(A,2);
                        end
                    end
                    
                    %%% Z-score
                    if contains(IdChannel{i},'Grab')
                        A = find(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(:,1) < 0 ...
                            & Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(:,3) < 1);
                    else
                        A = find(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(:,1) > 0 ...
                            & Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(:,3) < 1);
                    end
                    if ~isempty(A)
                        [~,idx] = max(abs(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A,1)));
                        A = A(idx);
                        if A > 1 && A < size(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.val(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A-1:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.time(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A-1:A+1,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.val(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A-1:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.time(o,:) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A-1:A+1,2);
                        elseif A > 1 && A >= size(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.val(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A-1:A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.time(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A-1:A,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.val(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A-1:A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.time(o,1:2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A-1:A,2);
                        elseif A == 1 && A < size(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.val(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.time(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A:A+1,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.val(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A:A+1,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.time(o,2:3) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A:A+1,2);
                        elseif A == 1 && A >= size(Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.val(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.AUC.time(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.AUC{o}(A,2);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.val(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A,1);
                            selectedMeasurements.LeverExtension.(IdChannel{i}).zscore.Peaks.time(o,2) = ...
                                Measurements.LeverExtension.(IdChannel{i}).zscore.Peaks{o}(A,2);
                        end
                    end
                end
            end
            %% %% Compute the AUC and the Peaks (lever press (dipper on))
%            
            idx_AUC1 = find(t_trials >= -1);
            idx_AUC1= idx_AUC1 (1);
            idx_AUC2= find(t_trials >= 3);
            idx_AUC2= idx_AUC2 (1);
            idx_AUC= idx_AUC1:idx_AUC2;
%             idx_AUC1 = find(t_trials == 0);
%             idx_AUC = idx_AUC1:(idx_AUC1+(idx_AUC2-1));
            clear dummie idx_AUC1 idx_AUC2
            
            for i = 1:length(IdChannel)
                Measurements.DipperOn.(IdChannel{i}).raw.AUC = cell(size(t_Dipper,1),1);
                Measurements.DipperOn.(IdChannel{i}).raw.Peaks = cell(size(t_Dipper,1),1);
                selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.val = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.time = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.val = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.time = ones(size(t_Dipper,1),3)*nan;
                
                Measurements.DipperOn.(IdChannel{i}).zscore.AUC = cell(size(t_Dipper,1),1);
                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks = cell(size(t_Dipper,1),1);
                selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.val = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.time = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.val = ones(size(t_Dipper,1),3)*nan;
                selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.time = ones(size(t_Dipper,1),3)*nan;
            end
            % Divide the trace in each sector above or bellow 0. Get the 
            % AUC and peak for each section
            for o = 1:size(t_Dipper,1)
                for i = 1:length(IdChannel)
                    %%% Raw
                    tmp = Stim_data.DipperOn.dFF.(IdChannel{i}).baseline_corrected(o,idx_AUC);
                    A = tmp >= 0;
                    B = find(A == 1);
                    if ~isempty(B)
                        C = diff(B);
                        D = find(C ~= 1);
                        zero_crossing = unique(sort([1 B(1) B(D) B(D+1) B(end) length(tmp)]));
                        tmp_AUC = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.val = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.idx = ones(length(zero_crossing)-1,1)*nan;
                        for oo = 1:length(zero_crossing)-1
                            tmp2 = tmp(zero_crossing(oo)+1:zero_crossing(oo+1));
                            tmp_AUC(oo) = trapz(tmp2);
                            [~,idx] = max(abs(tmp2));
                            tmp_peaks.val(oo) = tmp2(idx);
                            tmp_peaks.idx(oo) = zero_crossing(oo)+idx-1;
                        end
                        tmp2 = temp_t(zero_crossing);
                        Measurements.DipperOn.(IdChannel{i}).raw.AUC{o} = ...
                            [tmp_AUC diff(tmp2)' tmp2(1:end-1)'];
                        Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o} = ...
                            [tmp_peaks.val temp_t(tmp_peaks.idx)'];
                    else
                        [~,idx] = max(abs(tmp));
                        tmp_peaks.val(oo) = tmp(idx);
                        tmp_peaks.idx(oo) = idx;
                        Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o} = ...
                            [tmp_peaks.val(oo) temp_t(tmp_peaks.idx(oo))'];
                        Measurements.DipperOn.(IdChannel{i}).raw.AUC{o} = ...
                            [trapz(tmp) temp_t(length(idx_AUC)) 0];
                    end
                    
                    %%% Zscore
                    tmp = Stim_data.DipperOn.dFF.(IdChannel{i}).zscore_baseline_corrected(o,idx_AUC);
                    A = tmp >= 0;
                    B = find(A == 1);
                    if ~isempty(B)
                        C = diff(B);
                        D = find(C ~= 1);
                        zero_crossing = unique(sort([1 B(1) B(D) B(D+1) B(end) length(tmp)]));
                        tmp_AUC = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.val = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.idx = ones(length(zero_crossing)-1,1)*nan;
                        for oo = 1:length(zero_crossing)-1
                            tmp2 = tmp(zero_crossing(oo)+1:zero_crossing(oo+1));
                            tmp_AUC(oo) = trapz(tmp2);
                            [~,idx] = max(abs(tmp2));
                            tmp_peaks.val(oo) = tmp2(idx);
                            tmp_peaks.idx(oo) = zero_crossing(oo)+idx-1;
                        end
                        tmp2 = temp_t(zero_crossing);
                        Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o} = ...
                            [tmp_AUC diff(tmp2)' tmp2(1:end-1)'];
                        Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o} = ...
                            [tmp_peaks.val temp_t(tmp_peaks.idx)'];
                    else
                        [~,idx] = max(abs(tmp));
                        tmp_peaks.val(oo) = tmp(idx);
                        tmp_peaks.idx(oo) = idx;
                        Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o} = ...
                            [tmp_peaks.val(oo) temp_t(tmp_peaks.idx(oo))'];
                        Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o} = ...
                            [trapz(tmp) temp_t(length(idx_AUC)) 0];
                    end
                end
            end
            % Select the relevant sections: negative-positive-negative for 
            % DLight,positive-negative-positive for GaCh 
            for o = 1:size(t_Dipper,1)
                for i = 1:length(IdChannel)
                    %%% Raw
                    if contains(IdChannel{i},'Grab')
                        A = find(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(:,1) < 0 ...
                            & Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(:,3) < 1);
                    else
                        A = find(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(:,1) > 0 ...
                            & Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(:,3) < 1);
                    end
                    if ~isempty(A)
                        [~,idx] = max(abs(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A,1)));
                        A = A(idx);
                        if A > 1 && A < size(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.val(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A-1:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.time(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A-1:A+1,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.val(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A-1:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.time(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A-1:A+1,2);
                        elseif A > 1 && A >= size(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.val(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A-1:A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.time(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A-1:A,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.val(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A-1:A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.time(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A-1:A,2);
                        elseif A == 1 && A < size(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.val(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.time(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A:A+1,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.val(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.time(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A:A+1,2);
                        elseif A == 1 && A >= size(Measurements.DipperOn.(IdChannel{i}).raw.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.val(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.AUC.time(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.AUC{o}(A,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.val(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).raw.Peaks.time(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).raw.Peaks{o}(A,2);
                        end
                    end
                    
                    %%% Z-score
                    if contains(IdChannel{i},'Grab')
                        A = find(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(:,1) < 0 ...
                            & Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(:,3) < 1);
                    else
                        A = find(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(:,1) > 0 ...
                            & Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(:,3) < 1);
                    end
                    if ~isempty(A)
                        [~,idx] = max(abs(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A,1)));
                        A = A(idx);
                        if A > 1 && A < size(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.val(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A-1:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.time(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A-1:A+1,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.val(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A-1:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.time(o,:) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A-1:A+1,2);
                        elseif A > 1 && A >= size(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.val(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A-1:A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.time(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A-1:A,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.val(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A-1:A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.time(o,1:2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A-1:A,2);
                        elseif A == 1 && A < size(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.val(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.time(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A:A+1,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.val(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A:A+1,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.time(o,2:3) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A:A+1,2);
                        elseif A == 1 && A >= size(Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o},1)
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.val(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.AUC.time(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.AUC{o}(A,2);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.val(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A,1);
                            selectedMeasurements.DipperOn.(IdChannel{i}).zscore.Peaks.time(o,2) = ...
                                Measurements.DipperOn.(IdChannel{i}).zscore.Peaks{o}(A,2);
                        end
                    end
                end
            end
            %% Plot aligned results
            
            for cond = 1
                Int_scale = [];
                Color_scale = [];
                %%% Raw
                Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).baseline_corrected;
                Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).baseline_corrected;
                plot_aligned_results(Signal1,Signal2,t_trials,...
                    selectedMeasurements.LeverExtension.dLight.raw.Peaks.val(:,2),Int_scale,...
                    Color_scale,Condition{cond},IdChannel,PATH2SAVE,show_plot,save_plot)
                %%% Z-score
                Signal1 = Stim_data.(Condition{cond}).dFF.(IdChannel{1}).zscore_baseline_corrected;
                Signal2 = Stim_data.(Condition{cond}).dFF.(IdChannel{2}).zscore_baseline_corrected;
                plot_aligned_results(Signal1,Signal2,t_trials,...
                    selectedMeasurements.LeverExtension.dLight.zscore.Peaks.val(:,2),Int_scale,...
                    Color_scale,Condition{cond},IdChannel,PATH2SAVE,show_plot,save_plot)
            end
            %% Find AUC for GACh with lever press latency > 2s
            
            %ACh value
            AUC_val_mainDip= selectedMeasurements.LeverExtension.Grab.raw.AUC.val (:,2);
            AUC_val_mainDip_mean= mean(AUC_val_mainDip, 'omitnan');
            AUC_val_2sLat_GACh= AUC_val_mainDip([latency.leverExtention2pressMore2s], :);
            AUC_val_2sLat_GACh_mean= mean(AUC_val_2sLat_GACh, 'omitnan');
            
            %DA peak amplitude to lever extension
            DA_peak_amp= selectedMeasurements.LeverExtension.dLight.raw.Peaks.val (:,2);
            DA_peak_amp_mean = mean (DA_peak_amp, 'omitnan');
            DA_peak_amp_2sLat= DA_peak_amp([latency.leverExtention2pressMore2s], :);
            DA_peak_amp_2sLat_mean= mean(DA_peak_amp_2sLat, 'omitnan');
            
            %Rebound AUC
            AUC_rebound= selectedMeasurements.LeverExtension.Grab.raw.AUC.val (:,3);
            AUC_rebound_mean = mean(AUC_rebound, 'omitnan');
            AUC_rebound_2sLPL= AUC_rebound([latency.leverExtention2pressMore2s], :);
            AUC_rebound_2sLPL_mean= mean(AUC_rebound_2sLPL, 'omitnan');
            
            %DA AUC to lever extension
            DA_AUC= selectedMeasurements.LeverExtension.dLight.raw.AUC.val (:,2);
            DA_AUC_mean= mean (DA_AUC, 'omitnan');
            DA_AUC_2sLat= DA_AUC([latency.leverExtention2pressMore2s], :);
            DA_AUC_2sLat_mean= mean(DA_AUC_2sLat, 'omitnan');
            
            %DA peak amplitude to lever press
            DA_LP_peak_amp= selectedMeasurements.DipperOn.dLight.raw.Peaks.val(:,2);
            DA_LP_peak_amp_mean= mean(DA_LP_peak_amp, 'omitnan');
            DA_LP_peak_amp_2sLat= DA_LP_peak_amp([latency.leverExtention2pressMore2s], :);
            DA_LP_peak_amp_2sLat_mean= mean (DA_LP_peak_amp_2sLat, 'omitnan');
            
            %DA AUC to lever press
            DA_LP_AUC= selectedMeasurements.DipperOn.dLight.raw.AUC.val(:,2);
            DA_LP_AUC_mean= mean(DA_LP_AUC, 'omitnan');
            DA_LP_AUC_2sLat= DA_LP_AUC([latency.leverExtention2pressMore2s], :);
            DA_LP_AUC_2sLat_mean= mean(DA_LP_AUC_2sLat, 'omitnan');
            
            
            % ACh dip min
            ACh_dip_min= selectedMeasurements.LeverExtension.Grab.raw.Peaks.val (:, 2);
            ACh_dip_min_mean= mean(ACh_dip_min, 'omitnan');
            ACh_dip_2sLat= ACh_dip_min([latency.leverExtention2pressMore2s], :);
            ACh_dip_2sLat_mean= mean(ACh_dip_2sLat, 'omitnan');
            
            %time (pause length). Lever extension
            AUC_time_mainDip= selectedMeasurements.LeverExtension.Grab.raw.AUC.time (:,2);
            AUC_time_mainDip_mean= mean(AUC_time_mainDip, 'omitnan');
            AUC_time_2sLat_GACh= AUC_time_mainDip([latency.leverExtention2pressMore2s], :);
            AUC_time_2sLat_GACh_mean= mean(AUC_time_2sLat_GACh, 'omitnan');
            
            %time (pause length). Lever press
            AUC_press_time_mainDip= selectedMeasurements.DipperOn.Grab.raw.AUC.time (:,2);
            AUC_press_time_mainDip_mean= mean(AUC_press_time_mainDip, 'omitnan');
            AUC_press_time_2sLat_GACh= AUC_press_time_mainDip([latency.leverExtention2pressMore2s], :);
            AUC_press_time_2sLat_GACh_mean= mean(AUC_press_time_2sLat_GACh, 'omitnan');
            
            AUC_press_mainDip= selectedMeasurements.DipperOn.Grab.raw.AUC.val (:,2);
            AUC_press_mainDip_mean= mean(AUC_press_mainDip, 'omitnan');
            AUC_press_2sLat_mainDip= AUC_press_mainDip([latency.leverExtention2pressMore2s], :);
            AUC_press_2sLat_mainDip_mean= mean(AUC_press_2sLat_mainDip, 'omitnan');
            
            AUC_press_minDip= selectedMeasurements.DipperOn.Grab.raw.Peaks.val (:,2);
            AUC_press_minDip_mean= mean(AUC_press_minDip, 'omitnan');
            AUC_press_2sLat_minDip= AUC_press_minDip([latency.leverExtention2pressMore2s], :);
            AUC_press_2sLat_minDip_mean= mean(AUC_press_2sLat_minDip, 'omitnan');
            
            
            
%             %% Calculate the correlation between AUC and peaks of both sensors
            for o = 1:length(norm)
                % For the main AUC (peak of dLight and dip of Grab)
                [Trials_corr.(norm{o}).AUC.r,Trials_corr.(norm{o}).AUC.p,Trials_corr.(norm{o}).AUC.fitline] = corr_analysis_plot...
                    (selectedMeasurements.LeverExtension.dLight.(norm{o}).AUC.val(:,2),...
                    selectedMeasurements.LeverExtension.Grab.(norm{o}).AUC.val(:,2),corr_type,...
                    'AUC dLight','AUC Grab',[norm{o},' AUC correlation for trials'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' AUC correlation for trials']);
                
                [Trials_corr.(norm{o}).dLightAUC_latency2press.r,Trials_corr.(norm{o}).dLightAUC_latency2press.p,...
                    Trials_corr.(norm{o}).dLightAUC_latency2press.fitline] = corr_analysis_plot...
                    (latency.leverExtention2press,selectedMeasurements.LeverExtension.dLight.(norm{o}).AUC.val(:,2),...
                    corr_type,'Latency to press','AUC dLight',[norm{o},' Correlation dLight AUC to latency to press'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' Correlation dLight AUC to latency to press']);
                
                [Trials_corr.(norm{o}).GrabAUC_latency2press.r,Trials_corr.(norm{o}).GrabAUC_latency2press.p,...
                    Trials_corr.(norm{o}).GrabAUC_latency2press.fitline] = corr_analysis_plot...
                    (latency.leverExtention2press,selectedMeasurements.LeverExtension.Grab.(norm{o}).AUC.val(:,2),...
                    corr_type,'Latency to press','AUC Grab',[norm{o},' Correlation Grab AUC to latency to press'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' Correlation Grab AUC to latency to press']);
                
                [Trials_corr.(norm{o}).peak.r,Trials_corr.(norm{o}).peak.p,Trials_corr.(norm{o}).peak.fitline] = corr_analysis_plot...
                    (selectedMeasurements.LeverExtension.dLight.(norm{o}).Peaks.val(:,2),...
                    selectedMeasurements.LeverExtension.Grab.(norm{o}).Peaks.val(:,2),corr_type,...
                    'Peak dLight','Peak Grab',[norm{o},' Peaks correlation for trials'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' Peaks correlation for trials']);
                
                [Trials_corr.(norm{o}).dLightpeak_latency2press.r,Trials_corr.(norm{o}).dLightpeak_latency2press.p,...
                    Trials_corr.(norm{o}).dLightpeak_latency2press.fitline] = corr_analysis_plot...
                    (latency.leverExtention2press,selectedMeasurements.LeverExtension.dLight.(norm{o}).Peaks.val(:,2),...
                    corr_type,'Latency to press','Peak dLight',[norm{o},' Correlation dLight peak to latency to press'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' Correlation dLight peak to latency to press']);
                
                [Trials_corr.(norm{o}).Grabpeak_latency2press.r,Trials_corr.(norm{o}).Grabpeak_latency2press.p,...
                    Trials_corr.(norm{o}).Grabpeak_latency2press.fitline] = corr_analysis_plot...
                    (latency.leverExtention2press,selectedMeasurements.LeverExtension.Grab.(norm{o}).Peaks.val(:,2),...
                    corr_type,'Latency to press','Peak Grab',[norm{o},' Correlation Grab peak to latency to press'],show_plot,save_plot,...
                    PATH2SAVE,[norm{o},' Correlation Grab peak to latency to press']);
            end
            close all
            
            %% ITI Analysis
            % Identify the index of stimulation times to exclude them for the peak
            % selection
            if ~isempty(data2load)
                idxStim = unique(sort([reshape(Stim_data.LeverExtension.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.LeverExtension.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.LeverExtension.idx(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.DipperOn.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.DipperOn.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.DipperOn.idx(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.Reward.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.Reward.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.Reward.idx(:,idx_Init+1:end),2))]));
            else
                idxStim = unique(sort([reshape(Stim_data.LeverExtension.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.LeverExtension.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.LeverExtension.idx(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.DipperOn.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.DipperOn.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.DipperOn.idx(:,idx_Init+1:end),2))]));
            end
            dummie = 1:N;
            idx2include = setdiff(dummie,idxStim);
            clear dummie
            for i = 1:length(IdChannel)
                mean_dFF.(IdChannel{i}).raw = mean(double(dFF.(IdChannel{i}).filt(idx2include)));
                std_dFF.(IdChannel{i}).raw = std(double(dFF.(IdChannel{i}).filt(idx2include)));
                mean_dFF.(IdChannel{i}).zscore = mean(double(dFF.(IdChannel{i}).zscorefilt(idx2include)));
                std_dFF.(IdChannel{i}).zscore = std(double(dFF.(IdChannel{i}).zscorefilt(idx2include)));
            end
            %% Get overall correlation during ITI only
            %%% Using raw data (filt version)
            signal1 = dFF.dLight.filt(idx2include);
            signal2 = dFF.Grab.filt(idx2include);
            save_name = 'Raw ITI correlation with Grab lag';
            
            % Compute the correlation when lagging Grab agaist dLight
            [~,Ovrl_ITIcorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            clear signal1 signal2
            %% Get overall correlation during trials only with press latency > 2s (looking at cue induced (lever extension)) correlation and not movement
            %%% Using raw data (filt version)
            signal1 = dFF.dLight.filt(idx2include);
            signal2 = dFF.Grab.filt(idx2include);

            Stim_data_LeverExtension_2sLPL= Stim_data.LeverExtension.idx([latency.leverExtention2pressMore2s], :);
            if ~isempty(data2load)
                idxstim_2sLPL = unique(sort([reshape(Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),1,...
                    size(Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),1)*size...
                    (Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.DipperOn.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.DipperOn.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.DipperOn.idx(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.Reward.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.Reward.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.Reward.idx(:,idx_Init+1:end),2))]));
            else
                idxstim_2sLPL = unique(sort([reshape(Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),1,...
                    size(Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),1)*size...
                    (Stim_data_LeverExtension_2sLPL(:,idx_Init+1:end),2)),...
                    reshape(Stim_data.DipperOn.idx(:,idx_Init+1:end),1,...
                    size(Stim_data.DipperOn.idx(:,idx_Init+1:end),1)*size...
                    (Stim_data.DipperOn.idx(:,idx_Init+1:end),2))]));
            end
            
            signal1 = dFF.dLight.filt(idxstim_2sLPL);
            signal2 = dFF.Grab.filt(idxstim_2sLPL);
            save_name = 'Raw trial correlation with Grab lag';
           
            
            % Compute the correlation when lagging Grab agaist dLight
            [~,Ovrl_Trialcorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            clear signal1 signal2
%% Get correlation between DA and ACh for the lever extension for trials w/ press latency >2s

            Stim_data_LeverExtension_2sLPL= Stim_data.LeverExtension.idx([latency.leverExtention2pressMore2s], :);
            if ~isempty(data2load)
                idxstim_LeverExtension= unique(sort(reshape(Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),1,...
                    size(Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),1)*size...
                    (Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),2))));
            else
                idxstim_LeverExtension= unique(sort(reshape(Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),1,...
                    size(Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),1)*size...
                    (Stim_data_LeverExtension_2sLPL(:,t_trials >= -2 & t_trials <= 2),2))));
            end
            
            
            signal1 = dFF.dLight.filt(idxstim_LeverExtension);
            signal2 = dFF.Grab.filt(idxstim_LeverExtension);
            save_name = 'Raw trial lever press correlation with Grab lag';
           
            
            % Compute the correlation when lagging Grab agaist dLight
            [~,Ovrl_LeverExtcorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            clear signal1 signal2
            
%% Get correlation between DA and ACh for the lever press (dipper on) for trials w/ press latency >2s

            Stim_data_LeverPress_2sLPL= Stim_data.DipperOn.idx([latency.leverExtention2pressMore2s], :);
            if ~isempty(data2load)
                idxstim_LeverPress= unique(sort(reshape(Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),1,...
                    size(Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),1)*size...
                    (Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),2))));
            else
                idxstim_LeverPress= unique(sort(reshape(Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),1,...
                    size(Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),1)*size...
                    (Stim_data_LeverPress_2sLPL(:,t_trials >= -1 & t_trials <= 2),2))));
            end
            
            
            signal1 = dFF.dLight.filt(idxstim_LeverPress);
            signal2 = dFF.Grab.filt(idxstim_LeverPress);
            save_name = 'Raw trial lever press correlation with Grab lag';
           
            
            % Compute the correlation when lagging Grab agaist dLight
            [~,Ovrl_LeverPresscorr.raw] = ovrl_corr_calculation(signal1,signal2,...
                time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
            clear signal1 signal2
            
            %% Compute a moving threshoding for selecting events
            % Define parameters for moving threshold
            baseline_window = 60; %Window to calculate the threshold
            baseline_dt = 10; %Moving steps for the threshold
            
            temp_t = time_dwn - time_dwn(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= baseline_window);
            idx_baseline_window = dummie(1);
            dummie = 1:length(temp_t);
            dummie = dummie(temp_t >= baseline_dt);
            idx_baseline_dt = dummie(1);
            clear dummie temp_t
            
            % Index for the moving steps
            ix1 = 1:idx_baseline_dt:N-idx_baseline_window;
            ix2 = idx_baseline_window:idx_baseline_dt:N;
            ix = [ix1' ix2'];
            if ix(end,2) < N
               [ix] = [ix;[N-idx_baseline_window N]];
            end
            
            % Compute the moving threshold
            for channel = 1:length(IdChannel)
                moving_mean.(IdChannel{channel}).raw = ones(1,N)*nan;
                moving_std.(IdChannel{channel}).raw = ones(1,N)*nan;
                moving_mean.(IdChannel{channel}).zscore = ones(1,N)*nan;
                moving_std.(IdChannel{channel}).zscore = ones(1,N)*nan;
            end
            for o = 1:size(ix,1)
                idx = ix(o,1):ix(o,2);
                for channel = 1:length(IdChannel)
                    if o == 1
                        %Using the mean
                        %Using the median
                        moving_mean.(IdChannel{channel}).raw(idx) = nanmedian...
                            (dFF.(IdChannel{channel}).filt(idx));
                        moving_std.(IdChannel{channel}).raw(idx) = mad...
                            (dFF.(IdChannel{channel}).filt(idx)).*1.4826;
                        moving_mean.(IdChannel{channel}).zscore(idx) = nanmedian...
                            (dFF.(IdChannel{channel}).zscorefilt(idx));
                        moving_std.(IdChannel{channel}).zscore(idx) = mad...
                            (dFF.(IdChannel{channel}).zscorefilt(idx)).*1.4826;
                    else
                        %Using the mean
                        %Using the median
                        moving_mean.(IdChannel{channel}).raw(idx(end-idx_baseline_dt)+1:end) = nanmedian...
                            (dFF.(IdChannel{channel}).filt(idx));
                        moving_std.(IdChannel{channel}).raw(idx(end-idx_baseline_dt)+1:end) = mad...
                            (dFF.(IdChannel{channel}).filt(idx)).*1.4826;
                        moving_mean.(IdChannel{channel}).zscore(idx(end-idx_baseline_dt)+1:end) = nanmedian...
                            (dFF.(IdChannel{channel}).zscorefilt(idx));
                        moving_std.(IdChannel{channel}).zscore(idx(end-idx_baseline_dt)+1:end) = mad...
                            (dFF.(IdChannel{channel}).zscorefilt(idx)).*1.4826;
                    end
                end
            end
            
            %% Get Dlight events out of the stimulations
            for z = 1:length(norm)
                % Select epochs above thresholds
                dummie = 1:N;
                if z == 1
                    idx4peaks = dummie(dFF.dLight.filt>(moving_mean.dLight.(norm{z})+(moving_std.dLight.(norm{z})*THRESHOLD.dLight.peak)));
                else
                    idx4peaks = dummie(dFF.dLight.zscorefilt>(moving_mean.dLight.(norm{z})+(moving_std.dLight.(norm{z})*THRESHOLD.dLight.peak)));
                end
                A = diff(idx4peaks);
                B = find(A~=1);
                C = [0 B length(idx4peaks)];
                DA_peaks.(norm{z}).values = ones(1,length(C)-1)*nan;
                DA_peaks.(norm{z}).idx = ones(1,length(C)-1)*nan;
                for pks = 2:length(C)
                    if z == 1
                        tmp1 = dFF.dLight.filt(idx4peaks(C(pks-1)+1:C(pks)));
                    else
                        tmp1 = dFF.dLight.zscorefilt(idx4peaks(C(pks-1)+1:C(pks)));
                    end
                    [val,idx] = max(tmp1);
                    actual_idx = idx+idx4peaks(C(pks-1)+1)-1;
                    DA_peaks.(norm{z}).values(pks-1) = val;
                    DA_peaks.(norm{z}).idx(pks-1) = actual_idx;
                end
                [~,idx,~] = intersect(DA_peaks.(norm{z}).idx,idxStim);
                DA_peaks.(norm{z}).values(idx) = [];
                DA_peaks.(norm{z}).idx(idx) = [];
                N_DA_peaks.(norm{z}) = length(DA_peaks.(norm{z}).idx);
                clear dummie
                
                % Obtain the candidate events:
                temp_t = time_dwn - time_dwn(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= 5);
                idx_5s = dummie(1);
                
                for sensor = 1:length(IdChannel)
                    DA_events.(IdChannel{sensor}).(norm{z}).data = ones(N_DA_peaks.(norm{z}),(idx_5s*2)+1)*nan;
                    DA_events.(IdChannel{sensor}).(norm{z}).idx = ones(N_DA_peaks.(norm{z}),(idx_5s*2)+1)*nan;
                end
                
                for pks = 1:N_DA_peaks.(norm{z})
                    for sensor = 1:length(IdChannel)
                        if DA_peaks.(norm{z}).idx(pks)-idx_5s > 0 && DA_peaks.(norm{z}).idx(pks)+idx_5s <= N
                            if z == 1
                                DA_events.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).filt(DA_peaks.(norm{z}).idx(pks)-idx_5s:DA_peaks.(norm{z}).idx(pks)+idx_5s);
                            else
                                DA_events.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).zscorefilt(DA_peaks.(norm{z}).idx(pks)-idx_5s:DA_peaks.(norm{z}).idx(pks)+idx_5s);
                            end
                            DA_events.(IdChannel{sensor}).(norm{z}).idx(pks,:) = ...
                                    DA_peaks.(norm{z}).idx(pks)-idx_5s:DA_peaks.(norm{z}).idx(pks)+idx_5s;
                        end
                    end
                end
                t_events = time_dwn(1:(idx_5s*2)+1)-time_dwn(idx_5s+1); %Time for our selected epochs
                
                DA_peaks.(norm{z}).values(isnan(DA_events.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                DA_peaks.(norm{z}).idx(isnan(DA_events.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                for sensor = 1:length(IdChannel)
                    DA_events.(IdChannel{sensor}).(norm{z}).data...
                        (isnan(DA_events.(IdChannel{sensor}).(norm{z}).data(:,1)),:) = [];
                    DA_events.(IdChannel{sensor}).(norm{z}).idx...
                        (isnan(DA_events.(IdChannel{sensor}).(norm{z}).idx(:,1)),:) = [];
                end
                N_DA_peaks.(norm{z}) = length(DA_peaks.(norm{z}).idx);
                
                % Compute the AUC for the peaks and dips
                center = idx_5s+1;
                window = 2;
                idx_window = sum(t_events >= 0 & t_events <= window);
                temp_initPeak_idx = ones(1,N_DA_peaks.(norm{z}))*nan;
                DA_events.dLight.(norm{z}).auc.val = ones(1,N_DA_peaks.(norm{z}))*nan;
                DA_events.dLight.(norm{z}).auc.time = ones(1,N_DA_peaks.(norm{z}))*nan;
                DA_events.dLight.(norm{z}).peaks.val = ones(1,N_DA_peaks.(norm{z}))*nan;
                DA_events.dLight.(norm{z}).peaks.time = ones(1,N_DA_peaks.(norm{z}))*nan;
                for o = 1:N_DA_peaks.(norm{z})
                    trace = DA_events.dLight.(norm{z}).data(o,center-round(idx_window/2):center+round(idx_window/2));
                    curve = trace - mean(moving_mean.dLight.(norm{z})(DA_events.dLight.(norm{z}).idx(o,:)));
                    ctr = round(idx_window/2)+1;
                    DA_events.dLight.(norm{z}).peaks.val(o) = trace(ctr);
                    tmp1 =  curve(1:ctr-1);
                    [~,ix1] = findpeaks(tmp1*-1,1:ctr-1,'MinPeakHeight',-1);
                    if ~isempty(ix1)
                        ix1 = max(ix1);
                        zerocross = find(curve(ix1:ctr-1) <= 0);
                        if ~isempty(zerocross)
                            ix1 = zerocross(end) + ix1 - 1;
                        end
                    else
                        zerocross = find(tmp1 <= 0);
                        if ~isempty(zerocross)
                            ix1 = zerocross(end);
                        else
                            [~,ix1] = min(tmp1);
                            ix1 = ix1(end);
                        end
                    end
                    tmp2 =  curve(ctr+1:end);
                    [~,ix2] = findpeaks(tmp2*-1,ctr+1:length(curve),'MinPeakHeight',-1);
                    if ~isempty(ix2)
                        ix2 = min(ix2);
                        zerocross = find(curve(ctr+1:ix2) <= 0);
                        if ~isempty(zerocross)
                            ix2 = zerocross(1) + ctr;
                        end
                    else
                        zerocross = find(tmp2 <= 0);
                        if ~isempty(zerocross)
                            ix2 = zerocross(1) + ctr;
                        else
                            [~,ix2] = min(tmp2);
                            ix2 = ix2(1)+ctr;
                        end
                    end
                    temp_initPeak_idx(o) = (center-round(idx_window\2))+ix1;
                    DA_events.dLight.(norm{z}).peaks.time(o) = temp_t(ctr-ix1);
                    DA_events.dLight.(norm{z}).auc.val(o) = trapz(curve(ix1:ix2));
                    DA_events.dLight.(norm{z}).auc.time(o) = temp_t(length(ix1:ix2));
                end
                
                DA_events.Grab.(norm{z}).auc = cell(N_DA_peaks.(norm{z}),1);
                DA_events.Grab.(norm{z}).peaks = cell(N_DA_peaks.(norm{z}),1);
                selectedDA_events.dLight.(norm{z}).auc = DA_events.dLight.(norm{z}).auc;
                selectedDA_events.dLight.(norm{z}).peaks = DA_events.dLight.(norm{z}).peaks;
                selectedDA_events.Grab.(norm{z}).auc.val = ones(N_DA_peaks.(norm{z}),3)*nan;
                selectedDA_events.Grab.(norm{z}).auc.time = ones(N_DA_peaks.(norm{z}),3)*nan;
                selectedDA_events.Grab.(norm{z}).peaks.val = ones(N_DA_peaks.(norm{z}),3)*nan;
                selectedDA_events.Grab.(norm{z}).peaks.time = ones(N_DA_peaks.(norm{z}),3)*nan;
                for o = 1:N_DA_peaks.(norm{z})
                    trace = DA_events.Grab.(norm{z}).data(o,temp_initPeak_idx(o):temp_initPeak_idx(o)+idx_window);
                    %                 trace = DA_events.Grab.data(o,center-idx_window:center+idx_window);
                    curve = trace - mean(moving_mean.Grab.(norm{z})(DA_events.Grab.(norm{z}).idx(o,:)));
                    A = curve >= 0;
                    B = find(A == 1);
                    if ~isempty(B)
                        C = diff(B);
                        D = find(C ~= 1);
                        zero_crossing = unique(sort([1 B(1) B(D) B(D+1) B(end) length(curve)]));
                        tmp_AUC = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.val = ones(length(zero_crossing)-1,1)*nan;
                        tmp_peaks.idx = ones(length(zero_crossing)-1,1)*nan;
                        for oo = 1:length(zero_crossing)-1
                            tmp2 = curve(zero_crossing(oo)+1:zero_crossing(oo+1));
                            tmp_AUC(oo) = trapz(tmp2);
                            [~,idx] = max(abs(tmp2));
                            tmp_peaks.val(oo) = tmp2(idx);
                            tmp_peaks.idx(oo) = zero_crossing(oo)+idx-1;
                        end
                        tmp2 = temp_t(zero_crossing);
                        DA_events.Grab.(norm{z}).auc{o} = ...
                            [tmp_AUC diff(tmp2)' tmp2(1:end-1)'-DA_events.dLight.(norm{z}).peaks.time(o)];
                        DA_events.Grab.(norm{z}).peaks{o} = ...
                            [tmp_peaks.val temp_t(tmp_peaks.idx)'-DA_events.dLight.(norm{z}).peaks.time(o)];
                    else
                        [~,idx] = max(abs(curve));
                        tmp_peaks.val(oo) = curve(idx);
                        tmp_peaks.idx(oo) = idx;
                        DA_events.Grab.(norm{z}).peaks{o} = ...
                            [tmp_peaks.val(oo) temp_t(tmp_peaks.idx(oo))'-DA_events.dLight.(norm{z}).peaks.time(o)];
                        DA_events.Grab.(norm{z}).auc{o} = ...
                            [trapz(curve) temp_t(length(idx_AUC)) DA_events.dLight.(norm{z}).peaks.time(o)*-1];
                    end
                    %                 [~,ix1] = findpeaks(curve*-1,'MinPeakHeight',1);
                    %                 DA_events.Grab.auc(o) = trapz(curve);
                end
                for o = 1:N_DA_peaks.(norm{z})
                    A = find(DA_events.Grab.(norm{z}).auc{o}(:,1) < 0 ...
                        & DA_events.Grab.(norm{z}).auc{o}(:,3) < 1);
                    if ~isempty(A)
                        [~,idx] = max(abs(DA_events.Grab.(norm{z}).auc{o}(A,1)));
                        A = A(idx);
                        if A > 1 && A < size(DA_events.Grab.(norm{z}).auc{o},1)
                            selectedDA_events.Grab.(norm{z}).auc.val(o,:) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A-1:A+1,1);
                            selectedDA_events.Grab.(norm{z}).auc.time(o,:) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A-1:A+1,2);
                            selectedDA_events.Grab.(norm{z}).peaks.val(o,:) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A-1:A+1,1);
                            selectedDA_events.Grab.(norm{z}).peaks.time(o,:) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A-1:A+1,2);
                        elseif A > 1 && A >= size(DA_events.Grab.(norm{z}).auc{o},1)
                            selectedDA_events.Grab.(norm{z}).auc.val(o,1:2) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A-1:A,1);
                            selectedDA_events.Grab.(norm{z}).auc.time(o,1:2) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A-1:A,2);
                            selectedDA_events.Grab.(norm{z}).peaks.val(o,1:2) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A-1:A,1);
                            selectedDA_events.Grab.(norm{z}).peaks.time(o,1:2) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A-1:A,2);
                        elseif A == 1 && A < size(DA_events.Grab.(norm{z}).auc{o},1)
                            selectedDA_events.Grab.(norm{z}).auc.val(o,2:3) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A:A+1,1);
                            selectedDA_events.Grab.(norm{z}).auc.time(o,2:3) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A:A+1,2);
                            selectedDA_events.Grab.(norm{z}).peaks.val(o,2:3) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A:A+1,1);
                            selectedDA_events.Grab.(norm{z}).peaks.time(o,2:3) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A:A+1,2);
                        elseif A == 1 && A >= size(DA_events.Grab.(norm{z}).auc{o},1)
                            selectedDA_events.Grab.(norm{z}).auc.val(o,2) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A,1);
                            selectedDA_events.Grab.(norm{z}).auc.time(o,2) = ...
                                DA_events.Grab.(norm{z}).auc{o}(A,2);
                            selectedDA_events.Grab.(norm{z}).peaks.val(o,2) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A,1);
                            selectedDA_events.Grab.(norm{z}).peaks.time(o,2) = ...
                                DA_events.Grab.(norm{z}).peaks{o}(A,2);
                        end
                    end
                end
                
                % Compute the correlations:
                [DA_events_corr.(norm{z}).AUC.r,DA_events_corr.(norm{z}).AUC.p,DA_events_corr.(norm{z}).AUC.fitline] = corr_analysis_plot...
                    (selectedDA_events.dLight.(norm{z}).auc.val',...
                    selectedDA_events.Grab.(norm{z}).auc.val(:,2),corr_type,...
                    'AUC dLight','AUC Grab',[norm{z},' AUC correlation for ITI dLight peaks'],show_plot,save_plot,...
                    PATH2SAVE,[norm{z},' AUC correlation for ITI dLight peaks']);
                [DA_events_corr.(norm{z}).peaks.r,DA_events_corr.(norm{z}).peaks.p,DA_events_corr.(norm{z}).peaks.fitline] = corr_analysis_plot...
                    (selectedDA_events.dLight.(norm{z}).peaks.val',...
                    selectedDA_events.Grab.(norm{z}).peaks.val(:,2),corr_type,...
                    'Peak dLight','Peak Grab',[(norm{z}),' Peaks correlation for ITI dLight peaks'],show_plot,save_plot,...
                    PATH2SAVE,[(norm{z}),' Peaks correlation for ITI dLight peaks']);
%                 close all
                
                % Plot aligned results:
                n_epochs = []; % Number of epochs to randomly select to plot. If empty all are used.
                Int_scale = []; % Limits for the y_axis in the traces plots. If empty, automatic selection.
                Color_scale = []; % Limits for the color axis in the heatmap plots plots. If empty, automatic selection.
                IdEvent = 'DA peaks'; % Name of the events used to align
                
                plot_ITI_aligned_results(DA_events.(IdChannel{1}).(norm{z}).data,DA_events.(IdChannel{2}).(norm{z}).data,t_events,DA_peaks.(norm{z}).values,n_epochs,Int_scale,...
                    Color_scale,IdEvent,IdChannel,norm{z},PATH2SAVE,show_plot,save_plot)
%                 close all
            end

            %% Get GRAB dips out of the stimulations
            for z = 1:length(norm)
                % Select epochs bellow thresholds
                dummie = 1:N;
                if z == 1
                    idx4dips = dummie(dFF.Grab.filt<(moving_mean.Grab.(norm{z})-(moving_std.Grab.(norm{z})*THRESHOLD.Grab.dip)));
                else
                    idx4dips = dummie(dFF.Grab.zscorefilt<(moving_mean.Grab.(norm{z})-(moving_std.Grab.(norm{z})*THRESHOLD.Grab.dip)));
                end
                A = diff(idx4dips);
                B = find(A~=1);
                C = [0 B length(idx4dips)];
                ACh_dips.(norm{z}).values = ones(1,length(C)-1)*nan;
                ACh_dips.(norm{z}).idx = ones(1,length(C)-1)*nan;
                for pks = 2:length(C)
                    if z == 1
                        tmp1 = -1*dFF.Grab.filt(idx4dips(C(pks-1)+1:C(pks)));
                    else
                        tmp1 = -1*dFF.Grab.zscorefilt(idx4dips(C(pks-1)+1:C(pks)));
                    end
                    if length(tmp1) > 3
                        [val,idx] = findpeaks(tmp1,'NPeaks',1);
                        if ~isempty(idx)
                            actual_idx = idx+idx4dips(C(pks-1)+1)-1;
                        else
                            [val,idx] = min(tmp1);
                            actual_idx = idx+idx4dips(C(pks-1)+1)-1;
                        end
                    else
                        [val,idx] = min(tmp1);
                        actual_idx = idx+idx4dips(C(pks-1)+1)-1;
                    end
                    
                    if z == 1
                        ACh_dips.(norm{z}).values(pks-1) = dFF.Grab.filt(actual_idx);
                    else
                        ACh_dips.(norm{z}).values(pks-1) = dFF.Grab.zscorefilt(actual_idx);
                    end
                    ACh_dips.(norm{z}).idx(pks-1) = actual_idx;
                end
                [~,idx,~] = intersect(ACh_dips.(norm{z}).idx,idxStim);
                ACh_dips.(norm{z}).values(idx) = [];
                ACh_dips.(norm{z}).idx(idx) = [];
                N_ACh_dips.(norm{z}) = length(ACh_dips.(norm{z}).idx);
                clear dummie
                
                % Obtain candidate epochs
                for sensor = 1:length(IdChannel)
                    ACh_events.dips.(IdChannel{sensor}).(norm{z}).data = ones(N_ACh_dips.(norm{z}),(idx_5s*2)+1)*nan;
                    ACh_events.dips.(IdChannel{sensor}).(norm{z}).idx = ones(N_ACh_dips.(norm{z}),(idx_5s*2)+1)*nan;
                end
                
                for pks = 1:N_ACh_dips.(norm{z})
                    for sensor = 1:length(IdChannel)
                        if ACh_dips.(norm{z}).idx(pks)-idx_5s > 0 && ACh_dips.(norm{z}).idx(pks)+idx_5s <= N
                            if z == 1
                                ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).filt(ACh_dips.(norm{z}).idx(pks)-idx_5s:ACh_dips.(norm{z}).idx(pks)+idx_5s);
                            else
                                ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).zscorefilt(ACh_dips.(norm{z}).idx(pks)-idx_5s:ACh_dips.(norm{z}).idx(pks)+idx_5s);
                            end
                            ACh_events.(IdChannel{sensor}).(norm{z}).idx(pks,:) = ...
                                ACh_dips.(norm{z}).idx(pks)-idx_5s:ACh_dips.(norm{z}).idx(pks)+idx_5s;
                        end
                    end
                end
                
                ACh_dips.(norm{z}).values(isnan(ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                ACh_dips.(norm{z}).idx(isnan(ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                for sensor = 1:length(IdChannel)
                    ACh_events.dips.(IdChannel{sensor}).(norm{z}).data...
                        (isnan(ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(:,1)),:) = [];
                    ACh_events.dips.(IdChannel{sensor}).(norm{z}).idx...
                        (isnan(ACh_events.dips.(IdChannel{sensor}).(norm{z}).data(:,1)),:) = [];
                end
                N_ACh_dips.(norm{z}) = length(ACh_dips.(norm{z}).idx);
                
                % Plot aligned results:
                n_epochs = []; % Number of epochs to randomly select to plot. If empty all are used.
                Int_scale = []; % Limits for the y_axis in the traces plots. If empty, automatic selection.
                Color_scale = []; % Limits for the color axis in the heatmap plots plots. If empty, automatic selection.
                IdEvent = 'ACh dips'; % Name of the events used to align
                

                plot_ITI_aligned_results(ACh_events.dips.(IdChannel{1}).(norm{z}).data,...
                    ACh_events.dips.(IdChannel{2}).(norm{z}).data,t_events,ACh_dips.(norm{z}).values,...
                    n_epochs,Int_scale,Color_scale,IdEvent,IdChannel,norm{z},PATH2SAVE,show_plot,save_plot)
                close all
% 
            end
                                  
            %% Get GRAB peaks out of the stimulations
            for z = 1:length(norm)
                % Select epochs bellow thresholds
                dummie = 1:N;
                if z == 1
                    idx4peaks = dummie(dFF.Grab.filt>(moving_mean.Grab.(norm{z})+(moving_std.Grab.(norm{z})*THRESHOLD.Grab.peak)));
                else
                    idx4peaks = dummie(dFF.Grab.zscorefilt>(moving_mean.Grab.(norm{z})+(moving_std.Grab.(norm{z})*THRESHOLD.Grab.peak)));
                end
                A = diff(idx4peaks);
                B = find(A~=1);
                C = [0 B length(idx4peaks)];
                ACh_peaks.(norm{z}).values = ones(1,length(C)-1)*nan;
                ACh_peaks.(norm{z}).idx = ones(1,length(C)-1)*nan;
                for pks = 2:length(C)
                    if z == 1
                        tmp1 = dFF.Grab.filt(idx4peaks(C(pks-1)+1:C(pks)));
                    else
                        tmp1 = dFF.Grab.zscorefilt(idx4peaks(C(pks-1)+1:C(pks)));
                    end
                    if length(tmp1) > 3
                        [val,idx] = findpeaks(tmp1,'NPeaks',1);
                        if ~isempty(idx)
                            actual_idx = idx+idx4peaks(C(pks-1)+1)-1;
                        else
                            [val,idx] = max(tmp1);
                            actual_idx = idx+idx4peaks(C(pks-1)+1)-1;
                        end
                    else
                        [val,idx] = max(tmp1);
                        actual_idx = idx+idx4peaks(C(pks-1)+1)-1;
                    end

                    if z == 1
                        ACh_peaks.(norm{z}).values(pks-1) = dFF.Grab.filt(actual_idx);
                    else
                        ACh_peaks.(norm{z}).values(pks-1) = dFF.Grab.zscore(actual_idx);
                    end
                    ACh_peaks.(norm{z}).idx(pks-1) = actual_idx;
                end
                [~,idx,~] = intersect(ACh_peaks.(norm{z}).idx,idxStim);
                ACh_peaks.(norm{z}).values(idx) = [];
                ACh_peaks.(norm{z}).idx(idx) = [];
                N_ACh_peaks.(norm{z}) = length(ACh_peaks.(norm{z}).idx);
                clear dummie
                
                % Obtain candidate epochs
                for sensor = 1:length(IdChannel)
                    ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data = ones(N_ACh_peaks.(norm{z}),(idx_5s*2)+1)*nan;
                    ACh_events.peaks.(IdChannel{sensor}).(norm{z}).idx = ones(N_ACh_peaks.(norm{z}),(idx_5s*2)+1)*nan;
                end
                
                for pks = 1:N_ACh_peaks.(norm{z})
                    for sensor = 1:length(IdChannel)
                        if ACh_peaks.(norm{z}).idx(pks)-idx_5s > 0 && ACh_peaks.(norm{z}).idx(pks)+idx_5s <= N
                            if z == 1
                                ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).filt(ACh_peaks.(norm{z}).idx(pks)-idx_5s:ACh_peaks.(norm{z}).idx(pks)+idx_5s);
                            else
                                ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(pks,:) = ...
                                    dFF.(IdChannel{sensor}).zscorefilt(ACh_peaks.(norm{z}).idx(pks)-idx_5s:ACh_peaks.(norm{z}).idx(pks)+idx_5s);
                            end
                            ACh_events.peaks.(IdChannel{sensor}).(norm{z}).idx(pks,:) = ...
                                ACh_peaks.(norm{z}).idx(pks)-idx_5s:ACh_peaks.(norm{z}).idx(pks)+idx_5s;
                        end
                    end
                end
                
                ACh_peaks.(norm{z}).values(isnan(ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                ACh_peaks.(norm{z}).idx(isnan(ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(:,1))) = [];
                for sensor = 1:length(IdChannel)
                    ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data...
                        (isnan(ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(:,1)),:) = [];
                    ACh_events.peaks.(IdChannel{sensor}).(norm{z}).idx...
                        (isnan(ACh_events.peaks.(IdChannel{sensor}).(norm{z}).data(:,1)),:) = [];
                end
                N_ACh_peaks.(norm{z}) = length(ACh_peaks.(norm{z}).idx);
                
                % Plot aligned results:
                n_epochs = []; % Number of epochs to randomly select to plot. If empty all are used.
                Int_scale = []; % Limits for the y_axis in the traces plots. If empty, automatic selection.
                Color_scale = []; % Limits for the color axis in the heatmap plots plots. If empty, automatic selection.
                IdEvent = 'ACh peaks'; % Name of the events used to align
                
                plot_ITI_aligned_results(ACh_events.peaks.(IdChannel{1}).(norm{z}).data,...
                    ACh_events.peaks.(IdChannel{2}).(norm{z}).data,t_events,ACh_peaks.(norm{z}).values,...
                    n_epochs,Int_scale,Color_scale,IdEvent,IdChannel,norm{z},PATH2SAVE,show_plot,save_plot)
                close all
            end
            %% Save
            if done == 0 || overwrite == 1
                Params.MICEId = mice_list(m).name;
                Params.SESSIONId = sessions(s).name;
                Params.THRESHOLD = THRESHOLD;
                
                save([PATH2SAVE,'Session_analysis.mat'],'time_dwn','dFF','t_Lever','t_Dipper','t_LeverPress',...
                't_head_entry','t_reward','latency','lag2plot','Ovrl_corr','Ovrl_LeverPresscorr','Ovrl_LeverExtcorr','t_trials','Stim_data','Measurements',...
                'selectedMeasurements','moving_mean','moving_std','idx2include','mean_dFF',...
                'std_dFF','t_events','DA_peaks','DA_events','selectedDA_events','DA_events_corr','ACh_dips',...
                'ACh_peaks','ACh_events','Params', 'StimData_2sLat_GACh','StimData_2sLat_GACh_mean','StimData_2sLat_DA',...
                'StimData_2sLat_DA_mean','StimData_GACh','StimData_GACh_mean','StimData_DA','StimData_DA_mean',...
                'mean_LPlatency','StimData_LP_2sLat_GACh_mean','StimData_LP_2sLat_GACh','StimData_LP_2sLat_DA','StimData_LP_2sLat_DA_mean',...
                'AUC_val_mainDip','AUC_val_2sLat_GACh','AUC_val_2sLat_GACh_mean','AUC_time_mainDip',...
                'AUC_time_2sLat_GACh','AUC_time_2sLat_GACh_mean','AUC_val_mainDip_mean','AUC_time_mainDip_mean',...
                'ACh_dip_min_mean','ACh_dip_min','ACh_dip_2sLat','ACh_dip_2sLat_mean','AUC_press_mainDip','AUC_press_mainDip_mean',...
                'AUC_press_2sLat_mainDip','AUC_press_2sLat_mainDip_mean','AUC_press_minDip','AUC_press_minDip_mean','AUC_press_2sLat_minDip',...
                'AUC_press_2sLat_minDip_mean','AUC_press_time_mainDip','AUC_press_time_mainDip_mean','AUC_press_time_2sLat_GACh',...
                'AUC_press_time_2sLat_GACh_mean','DA_peak_amp','DA_peak_amp_mean','DA_peak_amp_2sLat','DA_peak_amp_2sLat_mean',...
                'DA_AUC','DA_AUC_mean','DA_AUC_2sLat','DA_AUC_2sLat_mean','DA_LP_peak_amp','DA_LP_peak_amp_mean','DA_LP_peak_amp_2sLat',...
                'DA_LP_peak_amp_2sLat_mean','DA_LP_AUC','DA_LP_AUC_mean','DA_LP_AUC_2sLat','DA_LP_AUC_2sLat_mean',...
                'AUC_rebound','AUC_rebound_mean','AUC_rebound_2sLPL','AUC_rebound_2sLPL_mean')
            end
        end
    end
end


