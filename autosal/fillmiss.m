function out=fillmiss(in,miss1,delmiss1,fill)

%function out=fillmiss(in,miss1,delmiss1,fill)
% in is (m x n) matrix with missing data indicated by miss1.
% Missing data are filled with fill.
% delmiss1 is tolerance (+-) miss1+-delmiss1.
% If miss1 is nan, then delmiss1 is irrelevant.

out=in;
[m,n]=size(in);
indicate=isnan(miss1);
		if indicate==0
% Fill no data with fill.
	[ii,jj]=find(in>=miss1-delmiss1&in<=miss1+delmiss1);
	nmiss=length(ii);
	for miss=1:nmiss
	out(ii(miss),jj(miss))=fill;
	end 
		else
	indnan=isnan(in);
	[ii,jj]=find(indnan==1);
	nmiss=length(ii);
	for miss=1:nmiss
	out(ii(miss),jj(miss))=fill;
	end 
		end 




	
