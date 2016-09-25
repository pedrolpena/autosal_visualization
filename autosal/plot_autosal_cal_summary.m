sdir=dir('salts');
if size(sdir)~=0
   cd salts
end;

output_file = 'process_output';
%%%%%
%+
% Open existing salts data base if one exists in this directory
%
clear autosal_salts
if exist ('autosal_salts_db.mat','file')
    load autosal_salts_db.mat
else
    sprintf('Warning: no autosal_salts_db.mat file in this directory \n \t %s \nExiting now. \n',pwd)
    return
end
%-
%%%%%

if exist ('calibations_db.mat','file')
    load calibations_db.mat
else
    calibration_data = get_calibrations_from_salts2 (autosal_salts);
end

output = make_salinity_file5 (autosal_salts);

%%%%%
%+
% Note this plots the duplicates
%
[output, duplicate_S] = average_and_remove_duplicates_aoml4 (output, [output_file '_salts.txt']);
%-
%%%%%

%%%%%
%+
% Plot calibration data time series
%
min_date = min(calibration_data(1).txt_date);
max_date = max(calibration_data(end).txt_date);
figure;
for II = 1:length(calibration_data)
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
legend('Good Start','Bad Start','Good End','Bad End','Location','southeast')
set(gca,'Box','off');
H=gca;
axesposition=get(gca,'Position');
ylimits = get(gca,'ylim');
hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100,24) sw_sals(ylimits(2)/100,24)],'Yaxislocation','right','xtick',[],'box','off');
axes(H)
set(H,'Box','on')
grid
title('Standard Water Values (x100) vs Date (left psu)') 
%print -depsc plot/calibrations_vs_date.ps



%%%%%
%+
% Get rid of weird large station numbers for axis limits
%
temp = vertcat(calibration_data.txt_station_id);
temp(find(temp>=500)) = [];
min_station = min(temp);
max_station = max(temp);
%
figure;
for II = 1: length(calibration_data)
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
    hold on;
    min_station = min(min(calibration_data(II).txt_station_id(I_start)), min_station);
    max_station = max(max(calibration_data(II).txt_station_id(I_start)), max_station);

    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'ok', 'Markerfacecolor', 'k','MarkerSize',10);
    if ~isempty(I_start)
        temp = vertcat(calibration_data.txt_station_id);
        temp(find(temp==999)) = [];
        min_station = min(temp);
        max_station = max(temp);
        %min_station = min(min(calibration_data(II).txt_station_id(I_start)), min_station);
        %max_station = max(max(calibration_data(II).txt_station_id(I_start)), max_station);
    end
    
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
    min_station = min(min(calibration_data(II).txt_station_id(I_start)), min_station);
    max_station = max(max(calibration_data(II).txt_station_id(I_start)), max_station);

    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==4);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
    if ~isempty(I_start)
        temp = vertcat(calibration_data.txt_station_id);
        temp(find(temp==999)) = [];
        min_station = min(temp);
        max_station = max(temp);
        %min_station = min(min(calibration_data(II).txt_station_id(I_start)), min_station);
        %max_station = max(max(calibration_data(II).txt_station_id(I_start)), max_station);
    end
    
    %
    % plot these dgood values again so that their symbols appear on the top
    % of the plot
    %
    I_start = find(calibration_data(II).txt_samp_nbr==1001 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id( I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
    I_start = find(calibration_data(II).txt_samp_nbr==1000 & calibration_data(II).txt_qc==2);
    plot (calibration_data(II).txt_station_id(I_start), 100*calibration_data(II).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
    set(gca,'xlim',[min_station max_station])
end
xlabel('Station Number')
set(gca,'xlim',[min_station-0.1*(max_station-min_station) max_station+0.1*(max_station-min_station)])
legend('Good Start','Bad Start','Good End','Bad End','Location','southeast')
set(gca,'Box','off');
H=gca;
axesposition=get(gca,'Position');
ylimits = get(gca,'ylim');
hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100/2,24) sw_sals(ylimits(2)/100/2,24)],'Yaxislocation','right','xtick',[],'box','off');
axes(H)
set(H,'Box','on')
grid
title('Standard Water Values (x100) vs Station (left psu)') 
%print -depsc plot/calibrations_vs_station_num.ps

cd ..
