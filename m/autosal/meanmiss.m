function m = meanmiss(x)
% MEANMISS - column matrix means with missing data
%  MEANMISS(X) returns the means of the columns of X 
%  as a row vector, where missing data in X are encoded as NaN's.
%  Returns the same result as MEAN if there are no NaNs in X.
%
%

% 
% VERSION: 
%	$Revision: 1.3 $
%	$Date: 2008/02/18 19:37:08 $
%	$Id: meanmiss.m,v 1.3 2008/02/18 19:37:08 duncombe Exp duncombe $
%
% CHANGELOG: 
%

notvalid = isnan(x);
x(notvalid) = zeros(size(x(notvalid)));
m = sum(x)./sum(1-notvalid);

