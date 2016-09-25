function output = make_salinity_file4 (salts)
%
% function output_salts = make_salinity_file (salts)
%
% Creates salinity text file output variable "output"
% that contains all the variables n column order for the Clivar Repeat
% Hydrography program (in a slightly changed order that they should be able
% to handle!).  POutting the columns in this should order should make them
% equivalently readable by Carlos's CTD Toolbox as well.
% Output = 
%   column 1 = output_sta  
%   column 2 = output_niskin
%   column 3 = output_bottle (sample bottle number)
%   column 4 = output_twiceratio 
%   column 5 = output_salinity 
%   column 6 = output_quality
%   column 7 = output_cast 
%   column 8 = salinity correction already applied 
%
% Version 1 was specific for Clivar A10 cruise.  Fixed many miscellaneous station problems.
% Version 2: Feb 2012
%   Uses only salts.txt* values, not newtxt* values.
%   Specific for Clivar A10 cruise.  Fixed many miscellaneous station
%   problems.
% Version 3: removed station specific editing 
% Version 4: Adds column for salinity correction that is added when
%   correcting autosal drift.  Note this needs to have several programs upstream 
%   and downstream modified. Here I check to make sure the structured
%   variable salts.salinity_correction exists.  If it does then it is
%   added.
%

clear output_sta
output_sta = salts(1).txt_station_id;
for II = 2: length(salts);
    output_sta = [output_sta' salts(II).txt_station_id']';
end;

%output_cast = ones(length(output_sta),1);
output_cast = vertcat(salts.txt_castnumber);

clear output_niskin
output_niskin = salts(1).txt_niskin_bottle;
for II = 2: length(salts);
    output_niskin = [output_niskin' salts(II).txt_niskin_bottle']';
end;

clear output_bottle
output_bottle = salts(1).txt_samp_nbr;
for II = 2: length(salts);
    output_bottle = [output_bottle' salts(II).txt_samp_nbr']';
end;

clear output_twiceratio
output_twiceratio = salts(1).txt_cond_ratio_correct;
for II = 2: length(salts);
    output_twiceratio = [output_twiceratio' salts(II).txt_cond_ratio_correct']';
end;

clear output_salinity
output_salinity = salts(1).txt_salinity_correct;
for II = 2: length(salts);
    output_salinity = [output_salinity' salts(II).txt_salinity_correct']';
end;

clear output_quality
output_quality = salts(1).txt_qc;
for II = 2: length(salts);
    output_quality = [output_quality' salts(II).txt_qc']';
end;
%output_quality = 2*ones(length(output_sta),1);

if ~isempty(salts(1).txt_salinity_correction)
    clear output_salinity_correction
    output_salinity_correction = salts(1).txt_qc;
    for II = 2: length(salts);
        output_salinity_correction = [output_salinity_correction' salts(II).txt_salinity_correction']';
    end;
    output = [output_sta output_niskin output_bottle ...
        output_twiceratio output_salinity output_quality output_cast output_salinity_correction];
else
    output = [output_sta output_niskin output_bottle ...
        output_twiceratio output_salinity output_quality output_cast];
end

II = find(output_sta == 1000);
output_sta (II,:) = [];
output (II,:) = [];

II = find(output_sta == 1001);
output_sta (II,:) = [];
output (II,:) = [];

II = find(output_sta == -9);
if ~isempty(II)
    output_sta (II,:) = [];
    output (II,:) = [];
end

II = find(output_sta == 888);
if ~isempty(II)
    output_sta (II,:) = [];
    output (II,:) = [];
end

%
% finally resort station numbers in order 0 to XXX.
% Sorts by station number
%
[junk,I] = sort(output(:,1),1);
junk=output(I,:);
output = junk;
return;








