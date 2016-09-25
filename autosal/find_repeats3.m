function index_dups = find_repeats3 (output)
%
% Function finds the same station and niskin bottle number in output and
% calls these "duplicate" samples.
%
% want to find duplicates of 
% station number and niskin number
%
% II should be (number_duplicate, 2) array of the indeces of all duplicates
%
% Version 1:  written for Clivar output files that had station, niskin, etc.
% Version 2: written for AOML files that have station, niskin, etc.
%   Output = 
%       column 1 = output_sta  
%       column 2 = output_niskin
%       column 3 = output_bottle (sample bottle number or run order number)
%       column 4 = output_twiceratio 
%       column 5 = output_salinity 
%       column 6 = output_quality
%       column 7 = output_cast 
%       column 8 = salinity correction
%
I_column_station = 1;
I_column_niskin = 2;
I_column_castnumber = 7;

num_duplicates = 0;
index_dups = [];
for II = 1:length(output)
    i_check_sta = output(II,I_column_station);
    i_check_niskin = output(II,I_column_niskin);
    i_check_castnumber = output(II,I_column_castnumber);
    JJ = find(output(:,I_column_station) == i_check_sta & output(:,I_column_niskin) == i_check_niskin & output(:,I_column_castnumber) == i_check_castnumber);
    if (length(JJ) == 2)
        % found one duplicate       
        num_duplicates = num_duplicates + 1;
        index_dups(:,num_duplicates) = JJ;
    elseif (length(JJ) > 2)
        % found more than one duplicate in the same Niskin bottle
        for KK = 2: length(JJ);
           num_duplicates = num_duplicates + 1 ;
           index_dups(:, num_duplicates) = [JJ(1) JJ(KK)]';
        end;
    end;
end;
if isempty(index_dups);
    fprintf(1,'\nNo Duplicates Found\n\n');
else
    index_dups=unique(index_dups','rows');
end
return;


    
    
