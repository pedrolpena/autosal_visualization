function calibration_data = get_calibrations_from_salts2 (salts)
%
% function:  
%       calibration_data = get_calibrations_from_salts2 (salts)
%       
% Version 2:  Feb 2012
%   Changed to use value salts.txt* instead of salts.newtxt*
%
for II = 1: length(salts)
    %
    % Figure out indeces of the standards
    % Note the first standard must be labeled "1000" and the closing standard for
    % each run must be labelled "1001".
    %
    standard1_idx=find(salts(II).sample_num==1000);
    standard2_idx=find(salts(II).sample_num==1001);
    calibration_data(II).sample_num = salts(II).sample_num([standard1_idx' standard2_idx']');
    calibration_data(II).sample_id_str = salts(II).sample_id_str([standard1_idx' standard2_idx']');
    calibration_data(II).uncorr_ratio = salts(II).uncorr_ratio ([standard1_idx' standard2_idx']');
    calibration_data(II).correction = salts(II).correction ([standard1_idx' standard2_idx']');
    calibration_data(II).adj_ratio = salts(II).adj_ratio ([standard1_idx' standard2_idx']');
    calibration_data(II).date_time = salts(II).date_time ([standard1_idx' standard2_idx']');
    calibration_data(II).dat_qc = salts(II).dat_qc ([standard1_idx' standard2_idx']');    
    
    standard1_idx=find(salts(II).txt_station_id==1000);
    standard2_idx=find(salts(II).txt_station_id==1001);
    calibration_data(II).txt_date = salts(II).txt_date ([standard1_idx' standard2_idx']');
    calibration_data(II).txt_cond_ratio = salts(II).txt_cond_ratio ([standard1_idx' standard2_idx']');
    calibration_data(II).txt_cond_ratio_correct = salts(II).txt_cond_ratio_correct ([standard1_idx' standard2_idx']');
    calibration_data(II).txt_qc = salts(II).txt_qc ([standard1_idx' standard2_idx']');
    calibration_data(II).txt_samp_nbr = salts(II).txt_samp_nbr ([standard1_idx' standard2_idx']');
    
    %idx_notstandard=find(salts(II).txt_station_id~=1000 & salts(II).txt_station_id~=1001);
    idx_notstandard=[max(standard1_idx)+1:1:min(standard2_idx)-1];
    start_sta = min(salts(II).txt_station_id (idx_notstandard));
    stop_sta = max(salts(II).txt_station_id (idx_notstandard));
    calibration_data(II).txt_station_id = [start_sta*ones(1,length(standard1_idx)) stop_sta*ones(1,length(standard2_idx))]';
         
    %standard1_idx=find(salts(II).txt_station_id==1000);
    %standard2_idx=find(salts(II).txt_station_id==1001);
    %calibration_data(II).txt_date = salts(II).txt_date ([standard1_idx' standard2_idx']');
    %calibration_data(II).txt_cond_ratio = salts(II).txt_cond_ratio ([standard1_idx' standard2_idx']');
    %calibration_data(II).txt_cond_ratio_correct = salts(II).txt_cond_ratio_correct ([standard1_idx' standard2_idx']');
    %calibration_data(II).txt_qc = salts(II).txt_qc ([standard1_idx' standard2_idx']');
    calibration_data(II).txt_niskin_bottle = salts(II).txt_niskin_bottle ([standard1_idx' standard2_idx']');
    
    %start_sta = salts(II).txt_station_id (standard1_idx(end)+1);
    %stop_sta = salts(II).txt_station_id (standard2_idx(1)-1);
    %calibration_data(II).txt_station_id = [start_sta*ones(1,length(standard1_idx)) stop_sta*ones(1,length(standard2_idx))]';
end

return

