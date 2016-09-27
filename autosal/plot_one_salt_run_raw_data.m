function plot_one_salt_run_raw_data (salts)

%
% Lets do some plotting for each sample
%  There are now length(unique_sample_str) number of "bottles"
%
% Plots out all reads (as long as there are less than or equal to six reads
% for each bottle.  If there are more, the median and averages still work,
% but all the values will not be plotted.
%

old_ylimits = [0 0];
for II = 1: length(salts.unique_sample_str);
    temp=salts.unique_sample_str(II);
    
    JJ_match_dat = find(strcmp(temp{1}, salts.sample_id_str));% replaced strmatch with this equivalent Pedro Pena 9.23.16
    JJ_match_raw = find(strcmp(temp{1}, salts.raw_sample_id_str));    
%     JJ_match_dat = strmatch(salts.unique_sample_str(II), salts.sample_id_str,'exact');
%     JJ_match_raw = strmatch(salts.unique_sample_str(II), salts.raw_sample_id_str,'exact');

    %figure;
    makefigexact4(7,4);
        wysiwyg;
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==1)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==1))), 'b*','markersize',10,'linewidth',2)
        hold on
 %
 % What you would like is for the text file to give you the median of all
 % reads (after first reducing to only 3 reads for standards.
 %
        h=text (mean(JJ_match_raw), median(100*salts.raw_uncorr_ratio(JJ_match_raw)), ['mean = ' num2str(mean(100*salts.raw_uncorr_ratio(JJ_match_raw)),7) ' +/- ' ...
            num2str(std(100*salts.raw_uncorr_ratio(JJ_match_raw)),7)  '    median = ' num2str(median(100*salts.raw_uncorr_ratio(JJ_match_raw)),7)]);
 
        set(h,'Horizontalalignment','center')
        title([upper(strrep(salts.file,'_',' ')) ' : ' strrep(char(salts.unique_sample_str(II)),'_',' ')])
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==2)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==2))), 'rx','markersize',10,'linewidth',2)
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==3)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==3))), 'g.','markersize',10,'linewidth',2)
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==4)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==4))), 'c+','markersize',10,'linewidth',2)
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==5)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==5))), 'ms','markersize',10,'linewidth',2)
        plot(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==6)), 100*salts.raw_uncorr_ratio(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==6))), 'yd','markersize',10,'linewidth',2)
              
        %
        % now plot the averaged data as found in the *.dat file
        % 
        %
        for KK = 1: length(JJ_match_dat)
            if (salts.reading_num(JJ_match_dat(KK)) == 0)
                % Can plot the average for all the samples
                plot([JJ_match_raw(1) JJ_match_raw(end)], 100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))],'k')
                plot([JJ_match_raw(1) JJ_match_raw(end)], 100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'k--' )
                plot([JJ_match_raw(1) JJ_match_raw(end)], 100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'k--' )
                plot([JJ_match_raw(1) JJ_match_raw(end)], 100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*0.00002,'k:','linewidth',2)
                plot([JJ_match_raw(1) JJ_match_raw(end)], 100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*0.00002,'k:','linewidth',2)

            else
                % match that up with the raw reading_num for the true ranges
                % need to find the range in question
                if (salts.reading_num(JJ_match_dat(KK)) == 1)
                    match_read = 1;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'b','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'b--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'b--' ,'linewidth',2)
                elseif (salts.reading_num(JJ_match_dat(KK)) == 2)
                    match_read = 2;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'r','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'r--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'r--' ,'linewidth',2)
                elseif (salts.reading_num(JJ_match_dat(KK)) == 3)
                    match_read = 3;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'g','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'g--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'g--' ,'linewidth',2)
                elseif (salts.reading_num(JJ_match_dat(KK)) == 4)
                    match_read = 4;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'c','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'c--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'c--' ,'linewidth',2)
                elseif (salts.reading_num(JJ_match_dat(KK)) == 5)
                    match_read = 5;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'm','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'm--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'm--' ,'linewidth',2)
                elseif (salts.reading_num(JJ_match_dat(KK)) == 6)
                    match_read = 5;
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))], 'y','linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] + 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'y--' ,'linewidth',2)
                    plot([min(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read))) max(JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==match_read)))], ...
                        100*[salts.uncorr_ratio(JJ_match_dat(KK)) salts.uncorr_ratio(JJ_match_dat(KK))] - 100*salts.uncorr_ratio_stnd_dev(JJ_match_dat(KK)), 'y--' ,'linewidth',2)
                end
            end
        end
        
        ylimits = get(gca,'ylim');
        set(gca,'ylim', [ylimits(1)-0.1*diff(ylimits) ylimits(2)+0.1*diff(ylimits)]);

        ylimits = get(gca,'ylim');
        xlimits = get(gca,'xlim');

        if (mean(old_ylimits) > mean(ylimits))
            h = text (xlimits(1)+0.1*diff(xlimits), ylimits(1)+0.1*diff(ylimits), [' Previous S range ' num2str(sw_sals(old_ylimits(1)/100,24),6) ' to ' num2str(sw_sals(old_ylimits(2)/100,24),6)]);
            set(h,'color',[1 0 0], 'fontweight','bold');
        elseif (mean(old_ylimits) < mean(ylimits))
            h = text (xlimits(1)+0.1*diff(xlimits), ylimits(1)+0.1*diff(ylimits), [' Previous S range ' num2str(sw_sals(old_ylimits(1)/100,24),6) ' to ' num2str(sw_sals(old_ylimits(2)/100,24),6)]);
            set(h,'color',[0 0 1], 'fontweight','bold');
        else
            h = text (xlimits(1)+0.1*diff(xlimits), ylimits(1)+0.1*diff(ylimits), [' Previous S range ' num2str(sw_sals(old_ylimits(1)/100,24),6) ' to ' num2str(sw_sals(old_ylimits(2)/100,24),6)]);
            set(h,'color',[0 0 0], 'fontweight','bold');
        end 
        
        set(h,'Horizontalalignment','left');
        old_ylimits = ylimits;

        grid on;
%
% plot secondary axes on the left side of the plot
%
        set(gca,'Box','off');
        axesposition=get(gca,'Position');
        ylimits = get(gca,'ylim');
        hNewAxes = axes('Position', axesposition, 'color', 'none', 'Ylim', [sw_sals(ylimits(1)/100,24) sw_sals(ylimits(2)/100,24)], 'Yaxislocation', 'right', 'xtick', [], 'box', 'off');
        %pause;
        waitforbuttonpress;
        close;
end;


return;
