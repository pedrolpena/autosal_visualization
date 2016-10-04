% process_autosal
%
% process WBTS salts from Brown
%
% must be called while matlab is in the directory where the raw autosal
% data is located.  This is typically in the CTD data subdirectory called
% salts.
%
%%%%%%
%
% USAGE NOTES:
%   1. Update some fields in the *.m file like output_file, standard_sea_water_batch and save
%	process_autosal file as another file name (e.g. process_autosal_a16n.m)
%   2. output_file: defaults to process_output.txt
%   3. standard_sea_water_batch: Make sure your standard seawater batch number is listed below.
%      	If not, add a new batch number where XXX is and a ration value where YYY appears and remove the comment
%   	"%" symbol in the first column of the line.
%   4. i_want_many_plots =1 for plots, =0 for less plots
%   5. correction_method = 1 = (default)
%
%%%%%
%
clear;
apath;
workingDir=pwd();
cfgDir=[workingDir,filesep,'cfg'];

%=========================================================================
%=========================================================================
% This piece loads the contents config.cfg into
% cell matrix 'params'
% The file is first read and copied to a temporary file without spaces that
% could potentially confuse textscan.

fileToRead=[cfgDir,filesep,'config.cfg'];
fileToWrite=[cfgDir,filesep,'config.cfg.tmp'];
fid=fopen(fileToRead,'r');
fidW=fopen(fileToWrite,'w');

tline = fgetl(fid);
while ischar(tline)
    
    if ~ischar(tline)
        break;
    end
    
    if ~isempty(tline)
        fprintf(fidW,[tline,'\r\n']);
    end
    tline = fgetl(fid);
end

fclose(fidW);
fidW=fopen(fileToWrite,'r');
params=textscan(fidW,'%s %s','Delimiter','=','CommentStyle','#');

fclose(fid);
fclose(fidW);
delete(fileToWrite);
%======================================================================
%======================================================================
rawDataDir=get_cruise_variable_value(params,'rawDataDir');
databaseDir=get_cruise_variable_value(params,'databaseDir');
plotsDir=get_cruise_variable_value(params,'plotsDir');
logDir=get_cruise_variable_value(params,'logDir');

if strcmp(rawDataDir,'default')
    rawDataDir=[workingDir,filesep,'raw_data'];
end
if strcmp(databaseDir,'default')
    databaseDir=[workingDir,filesep,'database'];
end
if strcmp(plotsDir,'default')
    plotsDir=[workingDir,filesep,'plots'];
end
if strcmp(logDir,'default')
    logDir=[workingDir,filesep,'logs'];
end



if ~exist(cfgDir,'dir')
    mkdir(cfgDir);
end

if ~exist(rawDataDir,'dir')
    mkdir(rawDataDir);
end

if ~exist(databaseDir,'dir')
    mkdir(databaseDir);
end

if ~exist(plotsDir,'dir')
    mkdir(plotsDir);
end

if ~exist(logDir,'dir')
    mkdir(logDir);
end



% sdir=dir('salts');
% if size(sdir)~=0
%     cd salts
% end;
% delete('*.mat');
% UPDATE SECTION
output_file = 'process_output';
%======================================================================
%The following will read in the IAPSO Standard Seawater
%from a text file named "standard_sewater_batch_values.txt" the format
%has to be batch number [single space] K15 value.
%Example  "158 .99970"
%======================================================================

% fileToRead='standard_sewater_batch_values.txt';
% fid=fopen([cfgDir,filesep,fileToRead],'r');
% cruiseVars=textscan(fid,'%s %s','Delimiter',',','CommentStyle','#');
% fclose(fid);


fileToRead=[cfgDir,filesep,'standard_sewater_batch_values.txt'];
fileToWrite=[cfgDir,filesep,'standard_sewater_batch_values.txt.tmp'];
fid=fopen(fileToRead,'r');
fidW=fopen(fileToWrite,'w');

tline = fgetl(fid);
while ischar(tline)
    
    if ~ischar(tline)
        break;
    end
    
    if ~isempty(tline)
        fprintf(fidW,[tline,'\r\n']);
    end
    tline = fgetl(fid);
end

fclose(fidW);
fidW=fopen(fileToWrite,'r');
cruiseVars=textscan(fidW,'%s %s','Delimiter',',','CommentStyle','#');

fclose(fid);
fclose(fidW);
delete(fileToWrite);
%======================================================================




%
% change the plots if you want
% yes=1 (plot out each reading); no=0 (only plot summary and duplicates)
i_want_many_plots = str2num(get_cruise_variable_value(params,'show_all_plots'));
correction_method = str2num(get_cruise_variable_value(params,'correction_method'));
%
% End update section
%
%%%%%%
pause_after_plot = str2num(get_cruise_variable_value(params,'pause_after_plot'));




%%%%%%
%
% Open existing salts data base if one exists in this directory
%
clear autosal_salts
if exist([databaseDir,filesep,'autosal_salts_db.mat'],'file');
    load([databaseDir,filesep,'autosal_salts_db.mat']);
end
%
%%%%%%




%%%%%%
%
% Get a list of all *.dat files
% These are assumed to all be autosal output files.
% Files can have *any* filename structure.
% Assumes that all *.dat files have the corresponding *.hrd and *.raw files
% with the same root filename.
%
flist = dir([rawDataDir,filesep,'*.dat']);

if length(flist) < 1
    sprintf('Warning: no *.dat files in this directory \n \t %s \nExiting now. \n',rawDataDir)
    return
end
%
%%%%%%




%%%%%%
%
% load and process only new *.dat files
%
% First find out which files do not appear in autosal_salts_db.mat already
%   This is kept as an index to the filename that is now (as true=1)
%
%
index_new_files = zeros(length(flist),1);
if exist('autosal_salts','var')
    for II = 1: length(flist)
        [pathstr,this_name,ext] = fileparts(flist(II).name);
        for JJ = 1: length(autosal_salts)
            index_new_files(II) = max(strcmpi(this_name,autosal_salts(JJ).file),index_new_files(II));
        end
    end
else
    index_new_files = zeros(length(flist),1);
end
%
%%%%%%




%%%%%
% Now read the new files
%
% The routine read_autosal_dat_raw_mb (called below)
% creates the autosal_salts structure that looks like:
% tmpsalts (later to become 'salts')
%    file - one file name with *.dat file containing:
%    sample_num         - for # of total reads including standards
%    sample_id_str      - for # of total reads including standards
%    bottle_label_str   - for # of total reads including standards
%    sample_type        - for # of total reads including standards
%    reading_num        - for # of total reads including standards
%    bath_t_str         - for # of total reads including standards
%    uncorr_ratio       - for # of total reads including standards
%    uncorr_ratio_stnd_dev - for # of total reads including standards
%    correction         - for # of total reads including standards
%    adj_ratio          - for # of total reads including standards
%    calc_s             - for # of total reads including standards
%    calc_s_stnddev     - for # of total reads including standards
%    date_time          - for # of total reads including standards
%    dat_qc             - for # of total reads including standards
% and these elements:  one for each of 10 values per single "read"
%    raw_sample_id_str
%    raw_reading_num
%    raw_uncorr_ratio
%    raw_adj_ratio
%    raw_calc_s
%    raw_qc
%    unique_sample_str
% These values are what was in the original file computed by the autosal
% software.
%
clear tmpsalts
KK_index_of_new_files = find(index_new_files == 0);
num_new_files = length(KK_index_of_new_files);
if (~isempty(KK_index_of_new_files))
    for II = 1: length(KK_index_of_new_files)
        flist(KK_index_of_new_files(II)).name
        this_file = flist(KK_index_of_new_files(II)).name;
        tmpsalts(II) = read_autosal_dat_raw_mb ([rawDataDir,filesep,flist(KK_index_of_new_files(II)).name]);
    end
else
    sprintf('No new *.dat files to process in this directory \n \t %s \nExiting now. \n',rawDataDir)
    return
end
%
%%%%%%




%%%%%%
%
% Do some plotting for new files only
%
if num_new_files > 0
    if i_want_many_plots
        for II = 1: num_new_files
            plot_one_salt_run_raw_data(tmpsalts(II),pause_after_plot);
        end
    end
end
%
%%%%%%





%%%%%%
%
%   make_one_btl_value_no_txt2
%   Adds the following fields to theResult
%           txt_station_id: [48x1 double]
%        txt_niskin_bottle: [48x1 double]
%             txt_samp_nbr: [48x1 double]
%            txt_tank_temp: [48x1 double]
%           txt_cond_ratio: [48x1 double]
%                 txt_date: [48x1 double]
%                   txt_qc;
%
if num_new_files > 0
    clear tmp
    for II = 1: num_new_files
        tmp(II) = make_one_btl_value_no_txt3 (tmpsalts(II));  % adds the averages to the file as *txt* elements
    end
    tmpsalts = tmp;
end
%
%%%%%%




%%%%%%
%
% Process New Files only
%
if num_new_files > 0
    %%%%%%
    %
    % correct_autosal_drift in matlab structure
    % (don't read in and write out new files)
    % This uses the txt values already present in autosal_salts structure
    %
    clear tmp
    for II = 1: num_new_files
        % 2.*standard_sea_water_batch(tmpsalts(II).batch_number)
        BatchNum=num2str(tmpsalts(II).batch_number);
        K15=str2num(get_cruise_variable_value(cruiseVars,BatchNum));
        tmp(II) = correct_autosal_drift2_matlab_structure_no_file3 (2.* K15, tmpsalts(II), 1, 0, 0,plotsDir,pause_after_plot);
        % Note:
        % correct_autosal_drift2_matlab_structure
        % adds the following fields to theResult
        %    txt_cond_ratio_correct
        %    txt_salinity_correct
        %close all;
    end
    tmpsalts = tmp;
    clear tmp;
    %
    %%%%%%
    
    
    %%%%%
    %+
    % Now add the new stations to the existing autosal_salts structure
    %
    if (exist('autosal_salts','var'))
        II_index_existing_salt_files = length(autosal_salts);
        % now append the new salt values to the end of the existing ones
        autosal_salts(end+1:end+num_new_files) = tmpsalts;
    else
        autosal_salts = tmpsalts;
    end
    %-
    %%%%%
end
%-
%%%%%


%%%%%
%+
% Create new output file
% With the corrected newtxt values create a salt file
%
output = make_salinity_file4 (autosal_salts);
[output, duplicate_S] = average_and_remove_duplicates_aoml4 (output, [output_file '_salts.txt'],plotsDir,databaseDir,logDir);
%-
%%%%%


%%%%%
%+
% Get all calibration data
%
calibration_data = get_calibrations_from_salts2 (autosal_salts);

save([databaseDir,filesep,'autosal_salts_db.mat'], 'autosal_salts','-v6');
save([databaseDir,filesep,'calibrations_db.mat'], 'calibration_data','-v6');

%%%%%%%%%%%%%%%%%%%%% Create an autosal run summary %%%%%%%%%%%%%%%%%%%%%%%

good_standard_qc=find(vertcat(calibration_data.txt_qc)==2);
%edit this if Autosals are changed out
%good_standard_qc=find(vertcat(calibration_data.txt_qc)==2 & vertcat(calibration_data.txt_station_id)<6);

aux_cond_ratio=vertcat(calibration_data.txt_cond_ratio);
aux_station=vertcat(calibration_data.txt_station_id);
aux_time=vertcat(calibration_data.date_time);

good_cond_ratio=aux_cond_ratio(good_standard_qc);
good_sal=sw_sals(good_cond_ratio/2,24);
good_station=aux_station(good_standard_qc);
good_time=aux_time(good_standard_qc);


[mc, bc, rc, smc, sbc] = do_mbari_fits2 (good_time, good_cond_ratio, 1);


cond_cruise_drift=mc(1)*(good_time(end)-good_time(1));

[ms, bs, rs, sms, sbs] = do_mbari_fits2 (good_time, good_sal, 1);

sal_cruise_drift=abs(ms(1))*(good_time(end)-good_time(1));

%plot(good_station,good_cond_ratio,'*')
%hold on;
%plot(good_station(1),good_cond_ratio(1),'g*')
%plot(good_station(end),good_cond_ratio(end),'r*')

cast_number=vertcat(autosal_salts.txt_castnumber);
sample_type=vertcat(autosal_salts.txt_sample_type);

good_sample=length(find(sample_type==1 & cast_number ~=8));
dup_s=length(duplicate_S);
submit_salts=good_sample-dup_s;
dup_perc=(dup_s/good_sample)*100;

nf=fopen([logDir,filesep,'autosal_summary.txt'],'w');
fprintf(nf, ['Cruise summary for autosal\n']);
fprintf(nf, ['Autosal cruise drift: ',num2str(cond_cruise_drift), '\n']);
fprintf(nf, ['Autosal cruise drift (salinity): ',num2str(sal_cruise_drift), '\n']);
fprintf(nf, ['Total Salts run = ', num2str(good_sample), '\n']);
fprintf(nf, ['Total dup_S run = ', num2str(dup_s), '\n']);
fprintf(nf, ['Total Salts submitted = ', num2str(submit_salts), '\n']);
fprintf(nf, ['Duplicate percentage = ', num2str(dup_perc)]);
fclose(nf);

clear mc bc rc smc sbc ms bs rs sms sbs cond_cruise_drift sal_cruise_drift
clear good_standard_qc aux_cond_ratio aux_station aux_time good_cond_ratio
clear good_sal good_station good_time cast_number sample_type good_sample
clear dup_s submit_salts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot calibration data time series

min_date = min(calibration_data(1).txt_date);
max_date = max(calibration_data(end).txt_date);
%figure;
makefigexact4(7,4);
for II = 1: length(calibration_data)
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_date(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
    min_date = min(min(calibration_data(II).txt_date( I_start)), min_date);
    max_date = max(max(calibration_data(II).txt_date( I_start)), max_date);
    hold on;
    
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_date(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'ok', 'Markerfacecolor', 'k','MarkerSize',10);
    if ~isempty(I_start)
        min_date = min(min(calibration_data(II).txt_date( I_start)), min_date);
        max_date = max(max(calibration_data(II).txt_date( I_start)), max_date);
    end
    
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_date( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
    min_date = min(min(calibration_data(II).txt_date( I_start)), min_date);
    max_date = max(max(calibration_data(II).txt_date( I_start)), max_date);
    
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_date( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
    if ~isempty(I_start)
        min_date = min(min(calibration_data(II).txt_date( I_start)), min_date);
        max_date = max(max(calibration_data(II).txt_date( I_start)), max_date);
    end
    %
    % plot these good values again so that their symbols appear on the top
    % of the plot
    %
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_date(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_date( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
end
datetick
set(gca,'xlim',[min_date-0.1*(max_date-min_date) max_date+0.1*(max_date-min_date)])
legend('Good Start','Bad Start','Good End','Bad End','Location','northoutside','orientation','horizontal')
set(gca,'Box','off');
axesposition=get(gca,'Position');
ylimits = get(gca,'ylim');
axis0=gca;
hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100/2,24) sw_sals(ylimits(2)/100/2,24)],'Yaxislocation','right','xtick',[],'box','off');
set(hNewAxes,'Position',axesposition);
set(axis0,'Position',axesposition*.9);
set(axis0,'Position',axesposition);
print('-depsc', [plotsDir,filesep,'calibrations_vs_date.eps']);
%print plot/calibrations_vs_date.fig
%savefig('plot/calibrations_vs_date.fig')

%min_station = min(calibration_data(1).txt_station_id);
%max_station = max(calibration_data(1).txt_station_id);

% Get rid of weird large station numbers for axis limits
%
temp = vertcat(calibration_data.txt_station_id);
temp(find(temp>=500)) = [];
min_station = min(temp);
max_station = max(temp);

%figure;
makefigexact4(7,4);
fig1Pos=get(gcf,'Position')
for II = 1: length(calibration_data)
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
    hold on;
    
    
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'ok', 'Markerfacecolor', 'k','MarkerSize',10);
    
    
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
    
    
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
    
    
    %
    % plot these dgood values again so that their symbols appear on the top
    % of the plot
    %
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
end
legend('Good Start','Bad Start','Good End','Bad End','Location','northoutside','orientation','horizontal')
xlabel('Station Number')
set(gca,'xlim',[min_station-0.1*(max_station-min_station) max_station+0.1*(max_station-min_station)])
set(gca,'Box','off');


axesposition=get(gca,'Position');
ylimits = get(gca,'ylim');
axis1=gca;
hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100/2,24) sw_sals(ylimits(2)/100/2,24)],'Yaxislocation','right','xtick',[],'box','off');


set(axis1,'Position',axesposition*.9);
set(axis1,'Position',axesposition);
set(hNewAxes,'Position',axesposition*.9);
set(hNewAxes,'Position',axesposition);
%set(gcf,'Position',fig1Pos*.9);
%set(gcf,'Position',fig1Pos);
print('-depsc', [plotsDir,filesep,'calibrations_vs_station.eps']);
%savefig('plot/calibrations_vs_station.fig')


