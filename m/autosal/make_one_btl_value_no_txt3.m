function theResult = make_one_btl_value_no_txt3 (salts)
%
%  There are now length(unique_sample_str) number of "bottles"
%  Multiple "reads" are possible for each bottle.
%
% Add the following fields to theResult
%           txt_station_id: [numberStationsx1 double]
%        txt_niskin_bottle: [numberStationsx1 double]
%             txt_samp_nbr: [numberStationsx1 double]
%            txt_tank_temp: [numberStationsx1 double]
%           txt_cond_ratio: [numberStationsx1 double]
%                 txt_date: [numberStationsx1 double]
%                   txt_qc;
%           txt_castnumber;
%
%
theResult = salts;
% Want station number and niskin number for each bottle value.
% Right now that is *not* in the files generated by OS/the autosal
% interface unless its input as part of the bottle number field by the
% autosal operator.
%
% First, take a guess at the station number from the filename
%
clear temp_stations
index_sta_number = strfind (salts.file,'_');
%temp = strrep (salts.file,'_',' ');
%stations = str2num(temp(index_sta_number(1):end));
if isempty(index_sta_number)
    temp =salts.file;
    stations = str2num(temp(end));
else
    temp = strrep (salts.file,'_',' ');
    stations = str2num(temp(index_sta_number(1):end));
end;
index_start = 0;
if length(stations) > 1
    index_max_sta = floor(length(salts.unique_sample_str)/(max(stations)-min(stations)+1));
    index_start = index_start+1;
    for JJ = stations(1):stations(end)
        temp_stations(index_start:index_start+index_max_sta-1) = JJ*ones(1,index_max_sta);
        index_start=index_start+index_max_sta;
    end
    if length(temp_stations) < length(salts.unique_sample_str)
        temp_stations (length(temp_stations)+1:length(salts.unique_sample_str)) = max(stations)*ones(length(salts.unique_sample_str)-length(temp_stations),1);
    end
    temp_stations = temp_stations(:);
    theResult.txt_station_id = temp_stations;
elseif length(stations)>0
    theResult.txt_station_id = stations*ones(length(salts.unique_sample_str),1);
else
    theResult.txt_station_id = zeros(length(salts.unique_sample_str),1);
end

%
% By default make the cast number = 1
theResult.txt_castnumber = ones(length(salts.unique_sample_str),1);


for II = 1: length(salts.unique_sample_str);
    temp=salts.unique_sample_str(II);
    
    JJ_match_dat = find(strcmp(temp{1}, salts.sample_id_str));% replaced strmatch with this equivalent Pedro Pena 9.23.16
    JJ_match_raw = find(strcmp(temp{1}, salts.raw_sample_id_str));
    %JJ_match_dat = strmatch(salts.unique_sample_str(II), salts.sample_id_str, 'exact');
    %JJ_match_raw = strmatch(salts.unique_sample_str(II), salts.raw_sample_id_str, 'exact');

    % Note:
    % JJ_match_raw = indeces of salts.raw* files that match each unique "bottle" 
    % JJ_match_raw(find(salts.raw_reading_num(JJ_match_raw)==1)) matches
    %   "read number" 1 from that bottle, etc
    %

    %
    % Fill in sample number, tank temp, date of sample 
    %
    theResult.txt_samp_nbr(II)  = salts.sample_num(JJ_match_dat(1));
    theResult.txt_tank_temp(II) = strread(char(strrep(salts.bath_t_str(JJ_match_dat(1)),'C','')),'%d');
    theResult.txt_date (II)= salts.date_time(JJ_match_dat(1));

    %
    % Put 1000 for all starting standards
    % Put 1001 for all ending standards as the station number
    %
    theResult.txt_sample_bottle_num(II) = NaN;
    if theResult.txt_samp_nbr(II) == 1000
        theResult.txt_station_id(II) = 1000;
        theResult.txt_niskin_bottle(II) = 1000;
    elseif theResult.txt_samp_nbr(II) == 1001
        theResult.txt_station_id(II) = 1001;
        theResult.txt_niskin_bottle(II) = 1001;
    else        
        temptext = deblank(char(salts.bottle_label_str(JJ_match_dat(1))));
        if ~isempty(temptext)
            % Then there is something the user has input in this string.
            % See if you can make sense out of what it is.
            bottle_id_as_number = str2num(temptext);
            if exist('bottle_id_as_number')
                if (length(bottle_id_as_number) >= 3)
                    theResult.txt_niskin_bottle(II) = bottle_id_as_number(2);
                    theResult.txt_station_id(II) = bottle_id_as_number(1);
                    theResult.txt_sample_bottle_num(II) = bottle_id_as_number(3);
                elseif (length(bottle_id_as_number) == 2)                          
                    theResult.txt_niskin_bottle(II) = bottle_id_as_number(2);
                    theResult.txt_station_id(II) = bottle_id_as_number(1);
                else
                    theResult.txt_niskin_bottle(II) = bottle_id_as_number(1);
                end
                if (length(bottle_id_as_number) >= 2)
                    % then double check that there is not a cast number in
                    % there as well
                    iblank = strfind(temptext,' ');
                    if (iblank > 4)
                        theResult.txt_castnumber(II) = mod(bottle_id_as_number(1),10);
                        theResult.txt_station_id(II) = floor(bottle_id_as_number(1)/10);
                    end
                end
            end
        else
            theResult.txt_niskin_bottle(II) = salts.sample_num(JJ_match_dat(1));
        end
    end
    
    % My Calculation:
    %            std(salts.raw_uncorr_ratio (JJ_match_raw)
    %            median(salts.raw_uncorr_ratio (JJ_match_raw)
    %
    % Note:  the standard *.txt files we've produced include 2*the conductivity
    % ratio, even though the autosal interface spits out the conductivity ratio
    % itself.
    %
    temp = 2*salts.raw_uncorr_ratio (JJ_match_raw);
    if (max(diff(temp)) > 6.e-5)
        theResult.txt_qc(II) = 2;
    else
        theResult.txt_qc(II) = 2;
    end
    temp(find(abs(temp-median(temp))>2.7*std(temp))) = [];
    
    %
    % For Debugging print out the values here:
    %    fprintf (1,[' ' num2str(theResult.txt_station_id(II)) ' \t ' num2str(theResult.txt_niskin_bottle(II)) ' \t ' num2str(max(diff(temp))) ' \t ' num2str(theResult.txt_qc(II)) ' \n'])

    theResult.txt_cond_ratio(II) = median(temp);
    theResult.txt_cond_ratio_std(II) = std(salts.raw_uncorr_ratio (JJ_match_raw));

    theResult.txt_bottle_label_str(II) = salts.bottle_label_str(JJ_match_dat(1));
    theResult.txt_sample_id_str(II) = salts.sample_id_str(JJ_match_dat(1));
    theResult.txt_sample_type(II) = salts.sample_type(JJ_match_dat(1));

end

theResult.txt_station_id = theResult.txt_station_id(:);
theResult.txt_samp_nbr=theResult.txt_samp_nbr(:);
theResult.txt_tank_temp = theResult.txt_tank_temp(:);
theResult.txt_date = theResult.txt_date(:);
theResult.txt_sample_bottle_num = theResult.txt_sample_bottle_num(:);
theResult.txt_niskin_bottle = theResult.txt_niskin_bottle(:);
theResult.txt_qc = theResult.txt_qc(:);
theResult.txt_cond_ratio=theResult.txt_cond_ratio(:);
theResult.txt_cond_ratio_std =theResult.txt_cond_ratio_std(:);
theResult.txt_bottle_label_str = theResult.txt_bottle_label_str(:);
theResult.txt_sample_id_str = theResult.txt_sample_id_str(:);
theResult.txt_sample_type = theResult.txt_sample_type(:);

return;
