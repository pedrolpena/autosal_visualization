function salts  = read_autosal_dat_raw_mb(this_file)
%function salts  = read_autosal_dat(this_file)
% read_autosal_dat - Read the *.dat files from the NOAA Ship  Ron Brown autosal system.
% 
% CTD Calibration toolbox
% 
% INPUT: 
%   fileName: path and name of the *.dat file to be translated.
%   Will also read the *.raw file of the same root extension in the same
%   directory.  The *.dat filename must be what is input.
%
% OUTPUT:
%   salts: sample structure containing the averaged autosal analysis.
%
% DESCRIPTION:	
%

%
% CHANGELOG: 
%   08-Jul-2004, Version 1.0
%		* Initial version.
%   02-Oct-2011 - rewritten to make life easier for reading these very very
%   standard files.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this_file = 'a10_008_009.dat';

[pathstr,this_name,ext] = fileparts(this_file);
tmp_salts.file = this_name;


% Read entire file into vectors of cell arrays or doubles depending on format string.
% File looks like this:
% a               b           c           d           e               f           g               h                   i               j           k           l               m                   n                           o
% SampleNumber	SampleID	BottleLabel	SampleType	ReadingNumber	DateTime	BathTemperature	UncorrectedRatio	UncorrectedRatioStandDev	Correction	AdjustedRatio	CalculatedSalinity	CalculatedSalinityStandDev	Comments
%
%a                 b               c   d   e   f                   g   h           i           j           k           l           m          
%1	Calibration#1	p152	0   1-Oct-2011 05:03:42	24C	0.999761	0.000002	0.000049	0.999810	0.000000	0.0000	
%1	Calibration#1		0	1	1-Oct-2011 05:03:42	24C	0.999761	0.000002	0.000049	0.999810	0.000000	0.0000	
%
%a   b               c   d   e   f                   g   h           i           j           k           l           m
%1	A10_008_009#1		1		1-Oct-2011 05:08:39	24C	0.993584	0.000003	-0.000025	0.993559	34.746640	0.0001	
%1	A10_008_009#1		1	1	1-Oct-2011 05:08:50	24C	0.993573	0.000003	-0.000025	0.993548	34.746190	0.0001	
%1	A10_008_009#1		1	2	1-Oct-2011 05:09:36	24C	0.993591	0.000003	-0.000025	0.993565	34.746880	0.0001	
%1	A10_008_009#1		1	3	1-Oct-2011 05:10:20	24C	0.993590	0.000002	-0.000025	0.993564	34.746840	0.0001	

[tmp_salts.sample_num, tmp_salts.sample_id_str, tmp_salts.bottle_label_str, tmp_salts.sample_type, tmp_salts.reading_num, date_time_str, tmp_salts.bath_t_str, ...
    tmp_salts.uncorr_ratio, tmp_salts.uncorr_ratio_stnd_dev, tmp_salts.correction, tmp_salts.adj_ratio, tmp_salts.calc_s, tmp_salts.calc_s_stnddev comment_field] ...
    = textread(this_file,'%d %s %s %d %d %s %s %f %f %f %f %f %f %s', 'headerlines',1,'delimiter','\t');

tmp_salts.date_time = datenum(date_time_str);
%
% Each bottle samples gets a quality flag of 2 to start
%
tmp_salts.dat_qc = 2 * ones(size(tmp_salts.date_time,1),1);

%
% Get rid of odd " marks that can occasionally appear in the output of the OS processing script
%
tmp_salts.sample_id_str = strrep(tmp_salts.sample_id_str,'"','');
tmp_salts.bath_t_str = strrep(tmp_salts.bath_t_str,'"','');

%
% replacing the standardization bottle number with 1000 or 1001
%
tmp_salts.sample_num (find(tmp_salts.sample_type==0)) = 1000;
std_idx=find(tmp_salts.sample_num==1000);
end_idx = find(diff(find(tmp_salts.sample_num==1000))>1)+1;
tmp_salts.sample_num(std_idx(end_idx:end))=1001;
end_idx = find(tmp_salts.sample_num==1001);


%
% Read bottle standard number
%
temp_txt=tmp_salts.bottle_label_str(find(tmp_salts.sample_type==0 & tmp_salts.reading_num==0));
if strncmp(temp_txt(1),'p',1)
    temp_txt=strrep(temp_txt,'p','');
else
    temp_txt=strrep(temp_txt,'P','');
end;
temp_batch_number = str2num(char(temp_txt));
tmp_salts.batch_number=temp_batch_number(1);

%
% Use the change below if you would like to plot all calibration values on
% the same plot in routine plot_one_salt_run_data
%
% Change sample ids for starting calibration to be:
% 'Calibration#1'
% Change sample ids for ending calibration to be:
% 'Calibration#2'
%
%tmp_salts.sample_id_str(std_idx) = {'Calibration#1'};
%tmp_salts.sample_id_str(end_idx) = {'Calibration#2'};

%
%
% Read equivalent raw data file as well.
% File should be in same directory etc as the dat file.
%
raw_file = [pathstr,filesep,this_name '.raw'];
[tmp_salts.raw_sample_id_str, tmp_salts.raw_reading_num, tmp_salts.raw_uncorr_ratio, tmp_salts.raw_adj_ratio, tmp_salts.raw_calc_s] = textread(raw_file,'%s %d %f %f %f ', 'headerlines',1,'delimiter','\t ');
std_idx = find(tmp_salts.sample_num==1000);
 
%
% Each bottle samples gets a quality flag of 2 to start
%
tmp_salts.raw_qc = 2 * ones(size(tmp_salts.raw_uncorr_ratio,1),1);

%
% Get rid of odd " marks that can occasionally appear in the output of the OS processing script
%
tmp_salts.raw_sample_id_str = strrep(tmp_salts.raw_sample_id_str,'"','');

%
% Change sample ids for starting calibration to be:
% 'Calibration#1'
% Change sample ids for ending calibration to be:
% 'Calibration#2'
%
%std_idx=strmatch('Calibration', tmp_salts.raw_sample_id_str)
%tmp_salts.raw_sample_id_str(std_idx(1:find(diff(std_idx)>1))) = {'Calibration#1'};
%tmp_salts.raw_sample_id_str(std_idx(find(diff(std_idx)>1)+1:end)) =  {'Calibration#2'};
%

% 
% get unique sample id strings
%
JJ = 1;
tmp_salts.unique_sample_str(JJ) = tmp_salts.sample_id_str(1);
for II = 1: length(tmp_salts.sample_id_str)
%     if isempty(strmatch(tmp_salts.sample_id_str(II), tmp_salts.unique_sample_str))
      if ~strcmp(tmp_salts.sample_id_str(II), tmp_salts.unique_sample_str)% Pero Pena 9.23.16
        JJ = JJ + 1;
        tmp_salts.unique_sample_str(JJ) = tmp_salts.sample_id_str(II);
    end
end
tmp_salts.unique_sample_str = tmp_salts.unique_sample_str(:);

salts=tmp_salts;

return;

