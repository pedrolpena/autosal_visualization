

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


station_number_to_plot=[];
while isempty(station_number_to_plot)
    station_number_to_plot = input('Enter station number you would like to plot. \n Enter a negative number to exit :  ');
end

while (station_number_to_plot >=0)
    % find station number
    for II = 1: length(autosal_salts)
        if ~isempty(find(autosal_salts(II).txt_station_id == station_number_to_plot))
            index_of_station_to_plot = II;
            break
        end
    end
    
    %figure
    makefigexact4(7,4);
        subplot(2,1,1)
            plot (autosal_salts(index_of_station_to_plot).txt_date, 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio, 'sk', 'Markerfacecolor', 'k','MarkerSize',5)
            hold on;
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio(I_start), 'ob', 'Markerfacecolor', 'b','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            datetick
            %legend('All','Good Start','Bad Start','Good End','Bad End','Location','best')
            xlabel('Date')
            ylabel('100*Conductivity Ratio')
            title('Uncorrected Conductivity Ratio/Salinity Values')
            set(gca,'Box','off');
            axesposition=get(gca,'Position');
            ylimits = get(gca,'ylim');
            H=gca;
            hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100,24) sw_sals(ylimits(2)/100,24)],'Yaxislocation','right','xtick',[],'box','off');        
            axes(H)
            set(gca,'Box','on');
            grid
        subplot(2,1,2)
            plot (autosal_salts(index_of_station_to_plot).txt_date, 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio_correct, 'sk', 'Markerfacecolor', 'k','MarkerSize',5)
            hold on;
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio_correct(I_start), 'ob', 'Markerfacecolor', 'b','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio_correct(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio_correct(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_date(I_start), 100/2*autosal_salts(index_of_station_to_plot).txt_cond_ratio_correct(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            datetick
            %legend('All','Good Start','Bad Start','Good End','Bad End','Location','best')
            xlabel('Date')
            ylabel('100*Conductivity Ratio')
            title('Corrected Conductivity Ratio/Salinity Values')
            set(gca,'Box','off');
            axesposition=get(gca,'Position');
            ylimits = get(gca,'ylim');
            H=gca;
            hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100,24) sw_sals(ylimits(2)/100,24)],'Yaxislocation','right','xtick',[],'box','off');
            axes(H)
            set(gca,'Box','on');
            grid

    %figure;
    makefigexact4(7,4);
        subplot(2,1,1)
            plot (autosal_salts(index_of_station_to_plot).txt_salinity_correct, 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction, 'sk', 'Markerfacecolor', 'k','MarkerSize',5)
            hold on;
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_salinity_correct(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'ob', 'Markerfacecolor', 'b','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_salinity_correct(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            plot (autosal_salts(index_of_station_to_plot).txt_salinity_correct(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            plot (autosal_salts(index_of_station_to_plot).txt_salinity_correct(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            xlabel('Salinity, psu')
            ylabel('Salinity Corrention X 1000, psu')
            grid
            title([ 'Run including station ' num2str(station_number_to_plot) '     File: ' autosal_salts(index_of_station_to_plot).file])
        subplot(2,1,2)
            plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle, 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction, 'sk', 'Markerfacecolor', 'k','MarkerSize',5)
            hold on;
            set(gca,'xlim',[0 25])
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            if ~isempty(I_start) 
                plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'ob', 'Markerfacecolor', 'b','MarkerSize',10)
            end
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1000 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            if ~isempty(I_start) 
                plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'og', 'Markerfacecolor', 'g','MarkerSize',10)
            end
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==4);
            if ~isempty(I_start) 
              plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vy', 'Markerfacecolor', 'y','MarkerSize',10)
            end
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            if ~isempty(I_start) 
                plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            end
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==1001 & autosal_salts(index_of_station_to_plot).txt_qc==2);
            if ~isempty(I_start) 
                plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            end
            I_start = find(autosal_salts(index_of_station_to_plot).txt_station_id==888);
            if ~isempty(I_start) 
                plot (autosal_salts(index_of_station_to_plot).txt_niskin_bottle(I_start), 1000*autosal_salts(index_of_station_to_plot).txt_salinity_correction(I_start), 'vr', 'Markerfacecolor', 'r','MarkerSize',10)
            end
            xlabel('Niskin Bottle')
            ylabel('Salinity Corrention X 1000, psu')
            grid
        
    station_number_to_plot=[];
    while isempty(station_number_to_plot)
        station_number_to_plot = input('Enter station number you would like to plot. \n Enter a negative number to exit :  ');
    end
end
    
