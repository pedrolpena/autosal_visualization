function theResult  = correct_autosal_drift2_matlab_structure_no_file3 (nominal_std, salts, use_method, input_slope_m, input_offset_b)
% CORRECT_AUTOSAL_DRIFT2_MATLAB_STRUCTURE - correct for autosal drift with pre/post standards 
% 
% CTD Calibration toolbox
% 
% INPUT: 
%   nominal_cond_standard_water: scalar representing the nominal value for
%		standard water conductivity ratio.  This number changes
% 		for each batch of standard seawater. 
%       e.g. nominal_std=2*0.99981;
%       Note:  it is assumed that the "nominal value is the same units as provided in
% 		the input file.
%   salts: structure of input data from the autosal interface
%   outfile: path and name of the output file which will contain the 
%		corrected averaged bottle sample salinity.
%       e.g. outfile = 'a10_014_015_out2.txt';
%
% OUTPUT:
%   theResult: the data structure of corrected values.  This is
%		essentially a duplicate of the information in the infile with an
%		additional varialbe that is the corrected value.
%
% DESCRIPTION:	
% 	This function is specific to the autosal configuration found on the
%	Ron Brown during the 2011 field season.  Modifications may be
%	necessary, but please rename accordingly.
%   This is just like correct_autosal_drift2.m except does not read in txt
%   files seperately.  Assumes txt files have already been read in and
%   placed in structure.
%
% SEE ALSO:
%		read_autosal_dat_raw		
%
% CHANGELOG: 
% Sept/2012 - wrote this
% Feb 2012 - check to make sure plot subdirectory exists.  If not, create
%               it.
% August 2013 - revised to handle samples outside the time range of
% standards (such as if trace metals were drawn).
% August 2013 - V3 changed to allow you to pick
%   use_method = 1 = standard least squares fit of standard data with QC
%              = 2 = uses mean of standard water values before/after run
%              = 3 = uses median of standard water values before/after run
%              = 4 = uses beginning standard water values (constant offset)
%              = 5 = uses ending standard water values (constant offset)
%              = 6 = uses slope/offset values specified in input
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
% Uncomment for debugging
% nominal_std = 2*0.99981;
% salts = tmpsalts(end);
% use_method = 1; (default =1)
% input_slope_m = 0;
% input_offset_b = 0;

% outfile = 'test.lis'

warning off MATLAB:divideByZero
error_for_standards = 0.00003;
error_for_drift     = 0.00006;


%
% check that plot subdirectory exists.  If not, create it.
if ~exist('plot','dir')
    mkdir plot
end

fprintf (1, ['\n\n\n Processing ' salts.file '\n\n\n']);


%
% Figure out indeces of the standards
% Note the first standard must be labeled "1000" and the closing standard for
% each run must be labelled "1001".
%
standard1_idx=find(salts.txt_station_id==1000);
standard2_idx=find(salts.txt_station_id==1001);
standard1_idx_all = standard1_idx;
standard2_idx_all = standard2_idx;


%
% Give all standards a "bad" flag" - reset flag after QC
%
salts.txt_qc(standard1_idx) = 4;
salts.txt_qc(standard2_idx) = 4;


%%%%%%%%%
%
% QC standard values
% 
% First get only the best 3 values
%
while length(standard1_idx) > 3
    [value index_max]=max(abs(salts.txt_cond_ratio(standard1_idx') - median(salts.txt_cond_ratio(standard1_idx'))));
    standard1_idx(index_max(1)) = [];
end

%
% Next check that there is not still one value that is an outlier
% 
% Note cond_ratio is actually 2*C(S,T,p)/C(35,15(ITS-68),0)
% A ratio difference of 0.0001 is like 0.002 psu
%                    of 0.00002 is like 0.0004 psu
%                    of 0.00006 is like 0.0012 psu
%
if max(salts.txt_cond_ratio(standard1_idx))-min(salts.txt_cond_ratio(standard1_idx)) > error_for_standards
    fprintf (1,['CorrectAutosalDrift2MatlabStructureNoFile:  \n    WARNING (' salts.file ') :  \n    Range of Initial Standard water values (larger than 0.00002) \n   = ' num2str(std(salts.txt_cond_ratio(standard1_idx')),13) ...
            '\n     Removing outlier and recompute.\n \n'])
    % One or more values must be really bad.  First find out if its the
    % first value by checking if the last two values are within range of
    % each other.
    if (length(standard1_idx) >= 3)
        if abs(salts.txt_cond_ratio(standard1_idx(3)) - salts.txt_cond_ratio(standard1_idx(2))) < error_for_standards
            % then use only the last two values
            standard1_idx(1)=[];
        else
            % get rid of the outlier whereever it is
            [value index_max] = max(abs(salts.txt_cond_ratio(standard1_idx') - median(salts.txt_cond_ratio(standard1_idx'))));
            standard1_idx(index_max(1)) = [];
        end
    end
end

%
% Now do the same thing for the final standards
%
while length(standard2_idx) > 3
    [value index_max] = max(abs(salts.txt_cond_ratio(standard2_idx') - median(salts.txt_cond_ratio(standard2_idx'))));
    standard2_idx(index_max(1)) = [];
end
%
% Next check that there is not still one value that is an outlier
%
if max(salts.txt_cond_ratio(standard2_idx))-min(salts.txt_cond_ratio(standard2_idx)) > error_for_standards
    fprintf (1,['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') :  \n Range of Final Standard water values is (larger than 0.00002) \n   = ' num2str(std(salts.txt_cond_ratio(standard1_idx')),13) ...
            '\n Removing outlier and recompute.\n \n'])
    % Then check the last two values
    if (length(standard2_idx) >= 3)
        if abs(salts.txt_cond_ratio(standard2_idx(3))- salts.txt_cond_ratio(standard2_idx(2))) < error_for_standards
            % then use only the last two values
            standard2_idx(1)=[];
        else
            % get rid of the outlier whereever it is
            [value index_max]=max(abs(salts.txt_cond_ratio(standard2_idx') - median(salts.txt_cond_ratio(standard2_idx'))))
            standard2_idx(index_max(1)) = [];
        end
    end
end



%
% Set standards that were actually used to "2" = good
%
salts.txt_qc(standard1_idx) = 2;
salts.txt_qc(standard2_idx) = 2;



%
% Make a variable containing only the standards
% x = the date/time of the standards
% y = the uncorrected standard ratios
%
x = salts.txt_date ([standard1_idx' standard2_idx']');
y = salts.txt_cond_ratio ([standard1_idx' standard2_idx']');
x_all = salts.txt_date ([standard1_idx_all' standard2_idx_all']');
y_all = salts.txt_cond_ratio ([standard1_idx_all' standard2_idx_all']');

%
% Make sure they are column vectors
%
x=x(:);
y = y(:);
x_all = x_all(:);
y_all = y_all(:);

%
% make sure there are no NaNs
%
II = find(isnan(y));
x(II)=[]; y(II) = [];

%
% Determine the standard water correction that needs to be applied.
%
btl_correction   = zeros(size(salts.txt_date,1), size(salts.txt_date,2));
btl_correction2  = zeros(size(salts.txt_date,1), size(salts.txt_date,2));

if isempty(standard1_idx)
    fpprintf (1, ['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') :  \n No starting standard found\n'])
    % Look for ending standard
    if isempty(standard2_idx)
        fpprintf (1, ['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') : \n No standards at all are found.  Cannot adjust values.\n'])
        btl_correction = zeros(length(salts.txt_cond_ratio),1);
    else
        fpprintf (1, ['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') : \n Only ending Standard is found.  Can only adjust mean offset for all bottles.\n'])
        btl_correction  = (mean  (salts.txt_cond_ratio(standard2_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio),1);
        btl_correction2 = (median(salts.txt_cond_ratio(standard2_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio),1);
        %
        % Do some very simple error checking:
        %
        if (std(salts.txt_cond_ratio(standard2_idx')) > 0.00001*2)
             fprintf (1,['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') : \n Ending Standard water values variance too large \n = ' num2str(std(salts.txt_cond_ratio(standard2_idx')),13) '\n'])
        end
    end
else
    if isempty(standard2_idx)
       fprintf (1, ['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') : \n Only beginning Standard is found.  Can only adjust mean offset for all bottles.\n'])
       btl_correction  = (mean  (salts.txt_cond_ratio(standard1_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio),1);
       btl_correction2 = (median(salts.txt_cond_ratio(standard1_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio),1);
        %
        % Do some very simple error checking:
        %
        if (std(salts.txt_cond_ratio(standard2_idx')) > 0.00001*2)
             fprintf (1,['CorrectAutosalDrift2MatlabStructureNoFile:  WARNING (' salts.file ') : \n Ending Standard water values variance too large \n = ' num2str(std(salts.txt_cond_ratio(standard2_idx')),13) '\n'])
        end
    else
        % Can do a linear offset correction as a function of time
        %
        % fit a linear polynomialerror_for_standards
        %
        % Many ways to do this e.g.:
        %   [coeff, S, tmp_best_fit_Rsq, tmp_best_fit_f, t] = molly_polyfit (x, y, 1, 1, infile);
        %   calculate_and_plot_linear_fit2 (x, y, ['standards for ' infile], 'date', ' autosal ratio', [infile '.ps'], 1)
        %
        %   Most robust method computes all possible varieties: do_mbari_fits2 (below)
        %       Note that the third and final argument is the unit number to write
        %       all the results to.  Unit=1 is the screen.  Use this unless you
        %       have already opened a logfile (otherwise there will be an error message).
        %
        fprintf(1,'\nUsing ALL beginning and ending standards with least squares fit:\n');
        [m, b, r, sm, sb] = do_mbari_fits2 (x, y, 1);
        % 
        % Also could make one median start and stop standard:
        fprintf(1,'\nUsing Median of all beginning and ending standards only:\n');
        x2 = [ mean(salts.txt_date(standard1_idx')) mean(salts.txt_date(standard2_idx'))];
        y2 = [ median(salts.txt_cond_ratio(standard1_idx')) median(salts.txt_cond_ratio(standard2_idx'))];
        [m2, b2, r2, sm2, sb2] = do_mbari_fits2 (x2, y2, 1);
        
        % Also could make one mean start and stop standard:
        fprintf(1,'\nUsing Mean of all beginning and ending standards only:\n');
        x3 = [ mean(salts.txt_date(standard1_idx')) mean(salts.txt_date(standard2_idx'))];
        y3 = [ mean(salts.txt_cond_ratio(standard1_idx')) mean(salts.txt_cond_ratio(standard2_idx'))];
        [m3, b3, r3, sm3, sb3] = do_mbari_fits2 (x3, y3, 1);
        

        % Use the starting standard only
        %fprintf(1,'\nUsing Median of all beginning standards only:\n');
        m4 = 0;
        b4 = median(salts.txt_cond_ratio(standard1_idx'));

        % Also could make one mean start and stop standard:
        %fprintf(1,'\nUsing Median of all ending standards only:\n');
        m5 = 0;
        b5 = median(salts.txt_cond_ratio(standard2_idx'));

        % fill slope and offset values just in case there are NANs.
        % Uses a zero slope and mean offset, b.
        m =fillmiss (m, NaN, 10, 0.);
        m2 =fillmiss (m2, NaN, 10, 0.);
        m3 =fillmiss (m3, NaN, 10, 0.);
        b  = fillmiss  (b, NaN, 10, mean(y));
        b2 = fillmiss  (b2, NaN, 10, mean(y2));
        b3 = fillmiss  (b3, NaN, 10, mean(y3));       
        
       
        %
        %  for the standards, the correction should give you
        %  y_orig - yfit = nominal_stnd
        %  where yfit = m*x + b
        %  hence y_new = y_orig - correction
        %  where correction = (yfit - nominal_stnd)
        %  y_new = y_orig - (yfit-nominal_stnd) = y_orig - yfit + nominal_stnd
        %
        % btl_correction = m(5)* date_cond + b(5) - nominal_std;
        % After looking at the first set of standards through station 36,
        % it looks best to use a median fit for all.
        
        if use_method == 1
            slope_to_use = m(5);
            offset_to_use = b(5);
        elseif use_method == 2
            slope_to_use = m2(5);
            offset_to_use = b2(5);
        elseif use_method == 3
            slope_to_use = m3(5);
            offset_to_use = b3(5);
        elseif use_method == 4
            slope_to_use = m4;
            offset_to_use = b4;
        elseif use_method == 5
            slope_to_use = m5;
            offset_to_use = b5;
        elseif use_method == 6
            slope_to_use = input_slope_m;
            offset_to_use = input_offset_b;
        else
            slope_to_use = m(5);
            offset_to_use = b(5);
        end
            
        
        %
        % First correct the entire run
        %      
        btl_correction = slope_to_use * salts.txt_date + offset_to_use - nominal_std;
        method_used = use_method;
        
            %
            % Do some very simple error checking:
            %
            if (std(salts.txt_cond_ratio(standard1_idx')) > 0.00001*2)
                fprintf (1,['\nCorrectAutosalDrift2MatlabStructureNoFile: WARNING (' salts.file ') :  \n ' ...
                 'Beginning Standard water values stnd dev is larger than 0.00002 (0.0008 psu) \n ' ...
                 '= ' num2str(std(salts.txt_cond_ratio(standard1_idx')),13) ...
                 '\n Consider removing outlier and recompute.\n '])
                %btl_correction = m2(5)* salts.txt_date + b2(5) - nominal_std;
            end
            if (std(salts.txt_cond_ratio(standard2_idx')) > 0.00001*2)
                fprintf (1,['\nCorrectAutosalDrift2MatlabStructureNoFile: WARNING (' salts.file ') :  \n ' ...
                 'Ending Standard water values stnd dev is larger than 0.00002 (0.0008 psu)\n ' ...
                 '= ' num2str(std(salts.txt_cond_ratio(standard1_idx')),13) ...
                 '\n   Consider removing outlier and recompute.\n '])
                %btl_correction = m2(5)* salts.txt_date + b2(5) - nominal_std;
            end
            if ((max(btl_correction) - min(btl_correction)) > 0.00006)
                fprintf (1,['\nCorrectAutosalDrift2MatlabStructureNoFile: WARNING (' salts.file ') :  \n ' ...
                'Standard water values changed more than 0.00006 (0.0024 psu) during run.  \n '  ... 
                'Check autosal peformance, room temperature, calibration values and scatter before and after run.\n'])
                if use_method == 7
                    fprintf (1,['     Using ending standard as constant correction.\n'])
                    btl_correction = (median  (salts.txt_cond_ratio(standard2_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio),1);
                    slope_to_use = 0;
                    offset_to_use = (median  (salts.txt_cond_ratio(standard2_idx')) - nominal_std);
                    method_used == 5;
                end
            end
            if ( abs(m2(5) - m(5)) > (abs(sm(5))) )
                fprintf (1,['\nCorrectAutosalDrift2MatlabStructureNoFile: WARNING (' salts.file ') :  \n     ' ...
                'Differences between slopes using Median and Mean reads is large.  \n     ' ...
                'Check calibration values scatter.  \n     ' ...
                'Consider removing outlier and recompute.\n ' ...
                '\n'])
                %btl_correction = m2(5)* salts.txt_date + b2(5) - nominal_std;
            end
 
        %
        % Then in case there are values outside the initial and final
        % standard (such as when analyzing trace metal salinities), use the
        % closest standard as a mean offset (no trend).
        LL = find( (salts.txt_date < min(salts.txt_date(standard1_idx_all'))) ); 
        if ~isempty(LL)
            btl_correction(LL)  = (median  (salts.txt_cond_ratio(standard1_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio(LL)),1);
        end
        MM = find( (salts.txt_date > max(salts.txt_date(standard2_idx_all'))) ); 
        if ~isempty(MM)
            btl_correction(MM)  = (median  (salts.txt_cond_ratio(standard2_idx')) - nominal_std) * ones(length(salts.txt_cond_ratio(MM)),1);
        end
        
        figure;
            subplot(2,1,1);
                plot(x_all,100*y_all,'ok','Markerfacecolor','k','MarkerSize',10)
                hold on;
                plot(x, 100*y,'ok','Markerfacecolor','b','MarkerSize',10)
                plot(x, 100*( m(5)*x+b(5)), 'b', 'linewidth', 2)
                plot(x, 100*(m2(5)*x+b2(5)), 'k--', 'linewidth', 1.5)
                plot(x, 100*(m3(5)*x+b3(5)), 'r:', 'linewidth', 1.5)
                plot(x, 100*(m3(5)*x+b3(5)), 'r:', 'linewidth', 1.5)
                plot(x, 100*(m4*x+b4), 'g:', 'linewidth', 1.5)
                plot(x, 100*(m5*x+b5), 'y:', 'linewidth', 1.5)
                legend ('All Standards', 'QC''d Standards', 'Fit from LSQ Reads', 'Fit from Median Reads', 'Fit from Mean Reads','Location','southeast') 
                datetick
                title (['Twice Standard Ratio times 100: ' upper(strrep(salts.file,'_',' '))])
  %              set(gca,'Box','off');
  %              axesposition=get(gca,'Position');
  %              ylimits = get(gca,'ylim');
  %              hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals(ylimits(1)/100/2,24) sw_sals(ylimits(2)/100/2,24)],'Yaxislocation','right','xtick',[],'box','off');

            subplot(2,1,2);
                plot(x_all, 100000*(y_all-(slope_to_use*x_all+offset_to_use)), 'ok','Markerfacecolor','k','MarkerSize',10)
                hold on;
                plot(x,     100000*(y              -(slope_to_use*x+offset_to_use)), 'ok','Markerfacecolor','b','MarkerSize',10)
                plot(x,     100000*((m(5)*x+b(5))  -(slope_to_use*x+offset_to_use)), 'b', 'linewidth', 2)
                plot(x,     100000*((m2(5)*x+b2(5))-(slope_to_use*x+offset_to_use)), 'k--', 'linewidth', 1.5)
                plot(x,     100000*((m3(5)*x+b3(5))-(slope_to_use*x+offset_to_use)), 'r--', 'linewidth', 1.5)
                plot(x,     100000*((m4*x+b4)-(slope_to_use*x+offset_to_use)), 'g--', 'linewidth', 1.5)
                plot(x,     100000*((m5*x+b5)-(slope_to_use*x+offset_to_use)), 'y--', 'linewidth', 1.5)
                datetick
                title ('Residual times 100000')
  %              set(gca,'Box','off');
  %              axesposition=get(gca,'Position');
  %              ylimits = get(gca,'ylim');
  %              hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals((nominal_std+ylimits(1)/100000)/2,24) sw_sals((nominal_std+ylimits(2)/100000)/2,24)],'Yaxislocation','right','xtick',[],'box','off');
                eval(['print -dpsc plot/cal_' salts.file '.ps'])
                %pause
                %waitforbuttonpress;
                
 %      figure;
 %           plot(salts.txt_date, 100*salts.txt_cond_ratio, 'ok','Markerfacecolor','b','MarkerSize',10)
 %           hold on;
 %           plot(salts.txt_date, 100*(salts.txt_cond_ratio-btl_correction),' ok','Markerfacecolor','r','MarkerSize',5)
 %           plot(x, 100*(y - (m(5)*x + b(5) - nominal_std)), ' *k','Markerfacecolor','g','MarkerSize',5) 
 %           datetick
 %           title (['Twice Standard Ratio times 100' upper(strrep(salts.file,'_',' '))])
 %           set(gca,'Box','off');
 %           axesposition=get(gca,'Position');
 %           ylimits = get(gca,'ylim');
 %           hNewAxes = axes('Position', axesposition, 'color','none','Ylim',[sw_sals((ylimits(1)/100)/2,24) sw_sals((ylimits(2)/100)/2,24)],'Yaxislocation','right','xtick',[],'box','off');
    end
end

btl_correction = btl_correction(:);
cond_ratio2_correct = salts.txt_cond_ratio-btl_correction;

output_correct=[salts.txt_station_id, salts.txt_niskin_bottle , ...
                salts.txt_samp_nbr  , cond_ratio2_correct     ,...
                sw_sals(cond_ratio2_correct/2, salts.txt_tank_temp)];

theResult = salts;
theResult.txt_cond_ratio_correct = cond_ratio2_correct;
theResult.txt_salinity_correct = sw_sals(cond_ratio2_correct/2,  24.);
theResult.txt_salinity_correction =  sw_sals(salts.txt_cond_ratio/2, 24.)  - theResult.txt_salinity_correct;
theResult.txt_salinity_method = method_used;

return

