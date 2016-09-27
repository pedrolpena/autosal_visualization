function the_returned_salts  = load_salts (cruiseid)
% load_salts : Load the autosal samples into the cruise database.
%
% function the_returned_salts  = load_salts(salts_file_name)
%
% CTD Calibration toolbox
% 
% INPUT: 
%   salts_file_name = character array containing the salt text file to
%   load.
%   cruiseid = character string that you 'registered' the cruise as (e.g.
%   ab1202)  If you don't know check getpref.m or register_cruise_ctd.m
%
% OUTPUT:
%   the_returned_salts: Combined salinity structure
%
%   Note: this also opens the cruise database (identified through the
%   cruiseid) and appends this salt structure to the database and saves it
%   on disk.
%
% DESCRIPTION:	
% This function reads the salinity samples from the matlab database.
%   1st column - Station Number; 
%   2nd - Niskin Bottle Number
%   3rd column - Sample bottle Number; 
%   4th- 2 * Conductivity ratio and 
%   5th column - Salinity
% 
% We assume that this data was already corrected for eventual drifts of
% the salinometer. The database is updated with the results.  Eventually it will become more cruise independent.
% Note that this can contain additional columns such as used by Repeat
% Hydrography, like cast number.  These columns are ignored here.
%
% Also assumes you have 
% 
% AUTHOR: 
%   Molly Baringer
%   Feb 2012
%
% CHANGELOG: 
%   08-Jul-2004, Version 1.0 Carlos Fonseca
%		* Initial version.
%   Feb 2012 - gutted this to make things much much easier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_cast_flag = 1;

sdir=dir('salts');
if size(sdir)~=0
   cd salts
end;
% Check for correct number of inputs.
error(nargchk(1,1,nargin))

% Retrieve the necessary preferences for this cruise
group = group_name(cruiseid);

% btl_dir is the calibration data directory 
if ispref(group,'cal_autosal_dir')
	cal_autosal_dir=getpref(group,'cal_autosal_dir');

    % cruise_dir is the calibration data directory 
    if ispref(group,'cal_data_dir')
        cruise_dir=getpref(group,'cal_data_dir');


        salt_db_name = (fullfile(cal_autosal_dir,'autosal_salts_db.mat'));

%
%       load the salts data structure
%
        if exist(salt_db_name,'file')
            eval(['load (''' salt_db_name ''');'])


            output = make_salinity_file4 (autosal_salts);
%           make_salinity_file2 (salts) creates one array with all standards removed:
%           Output = 
%               column 1 = output_sta  
%               column 2 = output_niskin
%               column 3 = output_bottle (sample bottle number)
%               column 4 = output_twiceratio 
%               column 5 = output_salinity 
%               column 6 = output_quality
%               column 7 = output_cast 
%
            [output, duplicate_S] = average_and_remove_duplicates_aoml4 (output, fullfile(cal_autosal_dir,'salts_as_loaded.txt'));
%           average_and_remove_duplicates averages all duplicate values and changes
%           the quality flag to duplicate = 5.
%           Writes out text file with the columns as listed above.
%

%
%           Now make sure the salts information is in the structure that the
%           CTDtoolbox is expecting
%
            the_returned_salts.station = output(:,1);
            the_returned_salts.niskin_bottle = output(:,2);
            the_returned_salts.sample_bottle = output(:,3);
            the_returned_salts.cond_ratio_2 = output(:,4);
            the_returned_salts.salinity = output(:,5);
            the_returned_salts.quality = output(:,6);
            the_returned_salts.cast = output(:,7);

            salts = the_returned_salts;
            
            % get rid of test salts
            if test_cast_flag
            test_sal=find(salts.station==0);
            salts.quality(test_sal)=9;
            end;
            
            close all
            
%           The cruise database.  Don't edit this line.
%
            db_file = fullfile(cruise_dir,[cruiseid '_db.mat']);

%
%           Open the cruise database and aapend or create this file with the salts
%           structure in it.
%
            if exist(db_file)==2
%                eval(['s=whos(''-file'',''' db_file ''');'])
%                II = strmatch('salts',char(s.name));
%                if ~isempty(II)
%                   load (db_file);
%                   eval(['save ' db_file ' ' sprintf('%s ',s.name)] );   
%                else
%                    eval(['save ' db_file ' ' sprintf('%s ',s.name) ' salts']);   
%                end
                save(db_file,'salts','-append');
            else
                save(db_file,'salts');
            end

            station_list = unique(vertcat(the_returned_salts.station));

            fprintf(1, 'Loaded salinity samples from station %f to %f \n', min(station_list), max(station_list));
            
        else
            sprintf('Warning no salinity database file \n\t%s  \n', salt_db)
            return
        end
    else
        errordlg(['The cal_data_dir preference was not set for cruise: ' cruiseid '.  Run register_cruise.m']);
        error(['The cal_data_dir preference was not set for cruise: ' cruiseid '.  Run register_cruise.m']);
    end
else
    errordlg(['The cal_autosal_dir preference was not set for cruise: ' cruiseid '.  Run register_cruise.m'])
    error(['The cal_autosal_dir preference was not set for cruise: ' cruiseid '.  Run register_cruise.m']);
end


if nargout==0
	assignin('caller','ans',salts);
    cd ..
    return
else
    varargout{1}=salts;
    cd ..
    return
end