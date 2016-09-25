
clear all;
pathn='./';
fileToRead='standard_sewater_batch_values.txt';
    fid=fopen([pathn fileToRead],'r');
    
    while ~feof(fid)
        linein=fgetl(fid);
        if ~isempty(linein)
            batch=textscan(linein,'%d %f','Delimiter',',','CommentStyle','%');
            standard_sea_water_batch(batch{1})=batch{2};
        end
    end

    fclose(fid);