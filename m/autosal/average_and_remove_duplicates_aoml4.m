function [theResult, duplicate_S] = average_and_remove_duplicates_aoml4 (output, outfile,plotsDir,databaseDir,logDir)
%
% Find duplicate samples in output array
%
% output = seven or eight column array of... 
%   [output_sta output_niskin output_bottle output_twiceratio output_salinity 
%          1         2            3             4                 5 
%                      output_quality output_cast (OPTIONAL) salinity_correction];
%                           6            7                         8
%
% this routine looks for duplicate output_niskin for each output_sta
% Averages these duplicates together provided they are less than 2*stddev
% of all duplicates away from the mean duplicate difference (i.e. tries to remove outliers
% before averaging).
%
% Output:
% outfile = text file of seven or eight columns as in output array, but with
% duplicates removed.
% Averaged duplicates qC flags are replaced with a "6"
%
% Version 1 called average_and_remove_duplicates for clivar with columns of
%   output different
% Version 2 called average_and_remove_duplicates_aoml for the format as
%   described above.
% Version 3 checks to see if output is 7 or 8 columns wide and deals with
%   the correct print statement for each case.
% Version4 only checks for eigth columns that also contain cast number
%   If cast number is present, then only averages duplicates if cast and
%   station number match. Saves the entire output row for bottle duplicates
%   so we can track down bad duplicates before they are averaged.
%
Icolumn_S = 5;
Icolumn_ratio = 4;
Icolumn_quality = 6;
Icolumn_salinity_correction = 8;
Icolumn_castnumber = 7;
Icolumn_station = 1;
Icolumn_niskin = 2;

[pathstr,this_name,ext] = fileparts(outfile);

duplicate_S = [];
index_dups = find_repeats3 (output);
% note that any index may appear more than once in this index_dups if more
% than one duplicate is drawn from a Niskin bottle

%%%%%
%+
% Store all the output data for every duplicate value
%
if ~isempty(index_dups)
saved_output_dups = output(sort([index_dups(:,1)' index_dups(:,2)']'),:);
save([databaseDir,filesep,'saved_output_dups.mat'],'saved_output_dups','index_dups','-v6');
end;
%-
%%%%%

if isempty(index_dups)
    fprintf(1,['\nNo Duplicates Found.  Writing ' outfile '\n\n']);
else
    duplicate_S     = zeros (size(index_dups,1), size(index_dups,2));
    duplicate_ratio = zeros (size(index_dups,1), size(index_dups,2));
    for II = 1:size(index_dups,1)
        duplicate_S (II,:) = output(index_dups(II,:),Icolumn_S);
        duplicate_ratio (II,:) = output(index_dups(II,:),Icolumn_ratio);
        duplicate_station (II,1)= output(index_dups(II),Icolumn_station);
    end
    if (size(output,2) == 8)
        duplicate_salinity_correction = zeros (size(index_dups,1), size(index_dups,2));
        for II = 1:size(index_dups,1)
            duplicate_salinity_correction(II,:) = output(index_dups(II,:),Icolumn_salinity_correction);
        end
    end
    if (size(output,2) == 9)
        duplicate_salinity_correction = zeros (size(index_dups,1), size(index_dups,2));
        for II = 1:size(index_dups,1)
            duplicate_castnumber(II,:) = output(index_dups(II,:),Icolumn_castnumber);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                CREATES A THE SALTS DUPLICATE TABLE                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1,
    k=1;
    staa=min(output(:,Icolumn_station));
    stab=max(output(:,Icolumn_station));
    for i=staa:1:stab,
        cast=unique(output(find(output(:,1)==i),Icolumn_castnumber));
        for c=1:length(cast),
            aux_sta=find(output(:,Icolumn_station)==i & output(:,Icolumn_castnumber)==cast(c));% & output(:,Icolumn_quality)==2);
            for j=1:24
                aux=find(output(aux_sta,Icolumn_niskin)==j);
                if length(aux)>1;
                    sta(k,1)=i;
                    niskin(k,1)=j;
                    botsal1(k,1)=output(aux_sta(aux(1)),5);
                    botsal2(k,1)=output(aux_sta(aux(2)),5);
                    k=k+1;
                end
            end
        end;
    end;
    fid=fopen([logDir,filesep,'duplicates_salt.tab'],'w');
    fprintf(fid,'%s','\begin{longtable}{ccccc}');
    fprintf(fid,'\n');    
    fprintf(fid,'%s','\caption{Duplicate salinity samples collected during the ABACO cruise.}\\');
    fprintf(fid,'\n');
    fprintf(fid,'%s','\hline');
    fprintf(fid,'\n');
    bp=[sprintf('%s\t','Station &','Niskin &','Salinity1 &','Salinity2 &','Differences \\')];
    fprintf(fid,bp);
    fprintf(fid,'\n');    
    fprintf(fid,'%s','\hline');
    fprintf(fid,'\n');   
    for ii=1:length(sta)
    fprintf(fid,'%2i\t%s%2i\t%s%8.3f\t%s%8.3f\t%s%8.3f%s\n',sta(ii), '&',...
        niskin(ii), '&', botsal1(ii), '&', botsal2(ii), '&', botsal1(ii)-botsal2(ii),'\\');
    end;
    fprintf(fid,'%s','\hline');
    fprintf(fid,'\n'); 
    fprintf(fid,'%s','\label{table:salt_dup}');
    fprintf(fid,'\n');          
    fprintf(fid,'%s','\end{longtable}');
    fclose(fid);
end;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    fonts=14;
    %
    % Average duplicates and plot differences
    %
    %    plot(mean(duplicate_S'),diff(duplicate_S'),'+')
    diff_dup_S = diff(duplicate_S');
    %Index_good_dups = find(abs(diff_dup_S) <= median(diff_dup_S)                  + 2*std(diff_dup_S));
    %Index_good_dups = find(abs(diff_dup_S) <= median(diff_dup_S(Index_good_dups)) + 2*std(diff_dup_S(Index_good_dups)));
    %Index_good_dups = find(abs(diff_dup_S) <= median(diff_dup_S(Index_good_dups)) + 2*std(diff_dup_S(Index_good_dups)));
    
    Index_good = find(diff_dup_S <= nanmedian(diff_dup_S) + 2*nanstd(diff_dup_S) & diff_dup_S >= nanmedian(diff_dup_S) - 2*nanstd(diff_dup_S));
    Index_good_dups = find(diff_dup_S <= nanmedian(diff_dup_S(Index_good)) + 2*nanstd(diff_dup_S(Index_good)) & diff_dup_S >= nanmedian(diff_dup_S(Index_good)) - 2*nanstd(diff_dup_S(Index_good)));

    %figure;
     makefigexact4(7,4);
     %makeaxis(2,2,6.5,9);
        number_dups = [1:length(diff_dup_S)];
        plot(number_dups, diff_dup_S,'+');
        hold on;
        plot(number_dups(Index_good_dups), diff_dup_S(Index_good_dups), 'r+');
        font(fonts)
        [xlimits]=get(gca,'xlim');
        plot(xlimits, [nanmedian(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))],'k--');
        plot(xlimits, [nanmedian(diff_dup_S(Index_good))+2*nanstd(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))+2*nanstd(diff_dup_S(Index_good))],'r:')
        plot(xlimits, [nanmedian(diff_dup_S(Index_good))-2*nanstd(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))-2*nanstd(diff_dup_S(Index_good))],'r:')
        [ylimits] = get(gca,'ylim');
        text([xlimits(1)+0.1*diff(xlimits)], [ylimits(2)-0.1*diff(ylimits)], ['Median = ' num2str(nanmedian(diff_dup_S(Index_good_dups)),5) ' +/- ' num2str(nanstd(diff_dup_S(Index_good_dups)),5)], 'fontsize', fonts)
        title (['Duplicate Salinity Differences from ' upper(strrep(this_name,'_',' '))], 'fontsize', 14);
        %set(gca,'xlim',[0.5 length(diff_dup_S)+0.5]);
        xlabel('Duplicate number');
        ylabel('Salinity, psu');
        set(gca,'xlim',[xlimits(1)-.03*xlimits(1),xlimits(2)+.03*xlimits(2)]);
        set(gca,'ylim',[ylimits(1)-.03*ylimits(1),ylimits(2)+.03*ylimits(2)]);
        %set(gca,'ylim',1.25*[ylimits]);

        print('-depsc', [plotsDir,filesep,'duplicates_salt_number.eps']);
        
     makefigexact4(7,4);
     %makeaxis(2,2,6.5,9);   
        plot(duplicate_station, diff_dup_S,'+');
        hold on;
        plot(duplicate_station(Index_good_dups), diff_dup_S(Index_good_dups), 'r+');
        font(fonts)
        [xlimits]=get(gca,'xlim');
        plot(xlimits, [nanmedian(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))],'k--');
        plot(xlimits, [nanmedian(diff_dup_S(Index_good))+2*nanstd(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))+2*nanstd(diff_dup_S(Index_good))],'r:')
        plot(xlimits, [nanmedian(diff_dup_S(Index_good))-2*nanstd(diff_dup_S(Index_good)) nanmedian(diff_dup_S(Index_good))-2*nanstd(diff_dup_S(Index_good))],'r:')
        [ylimits] = get(gca,'ylim');
        text([xlimits(1)+0.1*diff(xlimits)], [ylimits(2)-0.1*diff(ylimits)], ['Median = ' num2str(nanmedian(diff_dup_S(Index_good_dups)),5) ' +/- ' num2str(nanstd(diff_dup_S(Index_good_dups)),5)], 'fontsize', fonts)
        title ('Duplicate Salinity', 'fontsize', fonts);
        %set(gca,'ylim',[ylimits(1) ylimits(2)]);
        set(gca,'ylim',[ylimits(1)-.02*ylimits(1),ylimits(2)+.02*ylimits(2)]);
        %set(gca,'ylim',1.5*[ylimits])
        set(gca,'xlim',[xlimits(1)-.01*xlimits(1),xlimits(2)+.01*xlimits(2)]);         
        xlabel('Station number','fontsize', fonts);
        ylabel('Salinity Differences, psu','fontsize', fonts);
        print('-depsc', [plotsDir,filesep,'salts_duplicates.eps']);
    %
    % Average good duplicates and put the average value in the first
    % duplicate location and remove the second one
    %
    avg_S = meanmiss(duplicate_S')';
    avg_ratio = meanmiss(duplicate_ratio')';
    for II = 1: length(Index_good_dups);
        index_dups(Index_good_dups(II),1);
        JJ = find(index_dups(:,1) == index_dups(Index_good_dups(II),1));
        if (length(JJ) == 1)
            output (index_dups(Index_good_dups(II),1),Icolumn_ratio)   = avg_ratio(Index_good_dups(II));
            output (index_dups(Index_good_dups(II),1),Icolumn_S)       = avg_S(Index_good_dups(II));
            output (index_dups(Index_good_dups(II),1),Icolumn_quality) = 6;
            if (size(output,2) == 8)
                output (index_dups(Index_good_dups(II),1),Icolumn_salinity_correction) = meanmiss(duplicate_salinity_correction(Index_good_dups(II),:)')';
            end
        else
            % Then there is more than one row with duplicate values
            temp_index = unique(reshape(index_dups(JJ,:),size(JJ,1)*2,1));
            output (temp_index(1),Icolumn_ratio)   = meanmiss(output(temp_index,Icolumn_ratio));
            output (temp_index(1),Icolumn_S)       = meanmiss(output(temp_index,Icolumn_S));
            output (temp_index(1),Icolumn_quality) = 6;
        end
            
    end;
      
    %
    % Remove the duplicates (by convention remove the second duplicate).
    %
    output(index_dups(:,2),:) = [];
end


fid=fopen([logDir,filesep,outfile],'wt');
if (size(output,2) == 7)
    fprintf(fid,'%%Station\t Cast\t Niskin\t Sample_Bottle\t 2*Cond. Ratio Corrected Salinity\n');
elseif (size(output,2) == 8)
    fprintf(fid,'%%Station\t Cast\t Niskin\t Sample_Bottle\t 2*Cond. Ratio Corrected Salinity\n');
else
    fprintf(fid,'%%Station\t Cast\t Niskin\t Sample_Bottle\t 2*Cond. Ratio Corrected Salinity\n Salinity Correction \n');
end

for i = 1: size(output,1)
  fprintf (fid,'%11.5f', output(i,:));
  fprintf (fid,'\n');
end
      

fclose(fid);

fid=fopen([logDir,filesep,'cruise_salts.txt'],'wt');

%for a16s
for i = 1: size(output,1)
  fprintf(fid,'%%Station\t Cast\t Niskin\t Sample_Bottle\t 2*Cond. Ratio\t Corrected Salinity\t Salinity Quality\n');
  fprintf (fid,'%d\t%d\t%d\t%d\t%11.7f\t%11.7f\t%d', output(i,1),output(i,7),output(i,2),output(i,3),output(i,4),output(i,5),output(i,6));
  fprintf (fid,'\n');
end;

fclose(fid);

theResult = output;
return;


