function [m, b, r, sm, sb, rms_err, mean_per_err, mean_per_err2] = do_mbari_fits2 (volts, pegvtr, nf)
%
% function [m, b, r, sm, sb, rms_err, mean_per_err, mean_per_err2] = do_mbari_fits2 (volts, pegvtr, nf)
% Compute fits from Mbari
%   optionally writes to unit nf if nf .ne. 0
%
% Returns 1x7 arrays of fit values and statistics
%   1 = lsqfity = Model 1 least sq minimizing residuals in Y
%   2 = lsqfitxi = Model 1 least sq minimizing residuals in X
%   3 = lsqfityw = Model 1 least squares fit to WEIGHTED x,y-data pairs
%                 (used 10% of transport)
%   4 = lsqfityz = MINIMIZING the WEIGHTED residuals in Y only 
%                 (used 10% of transport)
%   5 = lsqfitma (Type II) = MINIMIZING the NORMAL deviates (passes thru
%   centroid)
%   6 = lsqfitgm (Type II) = GEOMETRIC MEAN
%       of the slopes from the regression of Y-on-X and X-on-Y
%   7 = lsqbisec (Type II) = LEAST SQUARES BISECTOR
%
% Output includes:
%     m = slope of fit to pegvtr = m*volts + b
%     b = intercept of fit to pegvtr = m*volts + b
%     r = correlation coefficient
%     sm = standard deviation of the slope error of fit to pegvtr = m*volts + b; m +/- sm
%     sb = standard deviation of the intercept error of fit to pegvtr = m*volts + b; b +/- sb
%
%     returned from compute_fit_stats
%         rms_err = root mean squared error (in Sv)
%         mean_per_err = % error relative to modelled or fit values, m*x+b
%         mean_per_err2 = % error relative to original series, y
%

m  = NaN * ones(1, 7);
b  = NaN * ones(1, 7);
r  = NaN * ones(1, 7);
sm = NaN * ones(1, 7);
sb = NaN * ones(1, 7);
rms_err = NaN * ones(1, 7);
mean_per_err = NaN * ones(1, 7);
mean_per_err2 = NaN * ones(1, 7);

%   [u,v,n,r,r2,YonLine,Ydiff] = model_II_reg(volts,pegvtr);
[my,by,ry,smy,sby]=lsqfity(volts,pegvtr);
fprintf(nf,['  lsqfity: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf);
m(1) = my;
b(1) = by;
r(1) = ry;
sm(1) = smy;
sb(1) = sby;
rms_err(1) = rmserr;
mean_per_err(1) = meanpererr;
mean_per_err2(1) = meanpererr2;

[my,by,ry,smy,sby]=lsqfitxi(volts,pegvtr);
fprintf(nf,['  lsqfitxi: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(2) = my;
b(2) = by;
r(2) = ry;
sm(2) = smy;
sb(2) = sby;
rms_err(2) = rmserr;
mean_per_err(2) = meanpererr;
mean_per_err2(2) = meanpererr2;

[my,by,ry,smy,sby]=lsqfityw(volts,pegvtr,0.1*meanmiss(pegvtr).*ones(size(pegvtr,1),size(pegvtr,2)));
fprintf(nf,['  lsqfityw err: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(3) = my;
b(3) = by;
r(3) = ry;
sm(3) = smy;
sb(3) = sby;
rms_err(3) = rmserr;
mean_per_err(3) = meanpererr;
mean_per_err2(3) = meanpererr2;

[my,by,ry,smy,sby]=lsqfityz(volts,pegvtr,0.1*meanmiss(pegvtr).*ones(size(pegvtr,1),size(pegvtr,2)));
fprintf(nf,['  lsqfityz: err \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(4) = my;
b(4) = by;
r(4) = ry;
sm(4) = smy;
sb(4) = sby;
rms_err(4) = rmserr;
mean_per_err(4) = meanpererr;
mean_per_err2(4) = meanpererr2;

fprintf(nf,['  ***** Type II \n']);

[my,by,ry,smy,sby]=lsqfitma(volts,pegvtr);
fprintf(nf,['  lsqfitma: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(5) = my;
b(5) = by;
r(5) = ry;
sm(5) = smy;
sb(5) = sby;
rms_err(5) = rmserr;
mean_per_err(5) = meanpererr;
mean_per_err2(5) = meanpererr2;

[my,by,ry,smy,sby]=lsqfitgm(volts,pegvtr);
fprintf(nf,['  lsqfitgm: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(6) = my;
b(6) = by;
r(6) = ry;
sm(6) = smy;
sb(6) = sby;
rms_err(6) = rmserr;
mean_per_err(6) = meanpererr;
mean_per_err2(6) = meanpererr2;

[my,by,ry,smy,sby]=lsqbisec(volts,pegvtr);
fprintf(nf,['  lsqbisec: \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;
m(7) = my;
b(7) = by;
r(7) = ry;
sm(7) = smy;
sb(7) = sby;
rms_err(7) = rmserr;
mean_per_err(7) = meanpererr;
mean_per_err2(7) = meanpererr2;

%[my,by,ry,smy,sby]=lsqcubic(volts,pegvtr,0.1*pegvtr.*ones(size(pegvtr,1),size(pegvtr,2)), 0.1*volts.*ones(size(pegvtr,1),size(pegvtr,2)),1e-5);
%fprintf(nf,['  lsqfityz: per err\t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
%[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;

%[my,by,ry,smy,sby]=lsqcubic(volts,pegvtr,0.1*meanmiss(pegvtr).*ones(size(pegvtr,1),size(pegvtr,2)), 0.1*meanmiss(volts).*ones(size(pegvtr,1),size(pegvtr,2)),1e-5);
%fprintf(nf,['  lsqfityz: err \t r=' num2str(ry) '\t r^2 = ' num2str(ry^2) '\t m = ' num2str(my) '+/-' num2str(smy) '\t b = ' num2str(by) '+/-' num2str(sby) '\n']);
%[rmserr, meanpererr, meanpererr2] = compute_fit_stats(volts, pegvtr, my, by, nf) ;

return;

function [rmserr, meanpererr, meanpererr2] = compute_fit_stats(X, Y, m, b, nf);

Ydiff = Y - b - m .* X;
rmserr = sqrt(sum(Ydiff .* Ydiff)/length(X));
meanpererr = meanmiss(100*abs(Ydiff)./ (m.*X + b));
meanpererr2 = 100*rmserr./meanmiss(Y);

if (nf ~= 0)  % 1=true 0=false
    fprintf(nf,['       new: \t rms=' num2str(rmserr) '\t %% error = ' num2str(meanpererr) ' RMS %% mean = ' num2str(meanpererr2) '\n']);
end

return

% lsqfity.m                                      by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
%
% M-file to calculate a "MODEL-1" least squares fit.
%
%     The line is fit by MINIMIZING the residuals in Y only.
%
%     The equation of the line is:     Y = my * X + by.
%
%     Equations are from Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       pp: 104, 108-109, 199.
%
%     Data are input and output as follows:
%
%         [my,by,ry,smy,sby] = lsqfity(X,Y)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%
%             my    =    slope
%             by    =    y-intercept
%             ry    =    correlation coefficient
%             smy   =    standard deviation of the slope
%             sby   =    standard deviation of the y-intercept

function [my,by,ry,smy,sby]=lsqfity(X,Y)

% Determine the size of the vector

n = length(X);

% Calculate the sums

Sx = sum(X);
Sy = sum(Y);
Sx2 = sum(X.^2);
Sxy = sum(X.*Y);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sx * Sy;
den = n * Sx2 - Sx^2;

% Calculate my, by, ry, s2, smy and sby

my = num / den;
by = (Sx2 * Sy - Sx * Sxy) / den;
ry = num / (sqrt(den) * sqrt(n * Sy2 - Sy^2));

diff = Y - by - my .* X;

s2 = sum(diff .* diff) / (n-2);
smy = sqrt(n * s2 / den);
sby = sqrt(Sx2 * s2 / den);

return;

% lsqfitxi.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2000 Jan 31.
%
% M-file to calculate a "MODEL-1" least squares fit.
%
%     The line is fit by MINIMIZING the residuals in X only.
%
%     The equation of the line is:     Y = mxi * X + bxi.
%
%     Equations are modified from those in Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       pp: 104, 108-109, 199.
%
%     Data are input and output as follows:
%
%         [mxi,bxi,rxi,smxi,sbxi] = lsqfitxi(X,Y)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%
%             mxi    =    slope
%             bxi    =    y-intercept
%             rxi    =    correlation coefficient
%             smxi   =    standard deviation of the slope
%             sbxi   =    standard deviation of the y-intercept

function [mxi,bxi,rxi,smxi,sbxi]=lsqfitxi(X,Y)

% Determine the size of the vector

n = length(X);

% Calculate the sums

Sx = sum(X);
Sy = sum(Y);
Sx2 = sum(X.^2);
Sxy = sum(X.*Y);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sy * Sx;
den = n * Sy2 - Sy^2;

% Calculate m, a, rx, s2, sm, and sb

mx = num / den;
a = (Sy2 * Sx - Sy * Sxy) / den;
rxi = num / (sqrt(den) * sqrt(n * Sx2 - Sx^2));

diff = X - a - mx .* Y;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sa = sqrt(Sy2 * s2 / den);

% Transpose coefficients

mxi = 1 / mx;
bxi = -a / mx;

smxi = mxi * sm / mx;
sbxi = abs(sa / mx);

return;


% lsqfityw.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
%
% M-file to calculate a "MODEL-1" least squares fit to WEIGHTED x,y-data pairs:
%
%     The line is fit by MINIMIZING the WEIGHTED residuals in Y only.
%
%     The equation of the line is:     Y = mw * X + bw.
%
%     Equations are from Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       for mw, bw, smw and sbw, see p. 98, example calculation in Table 6.2;
%       for rw, see p. 199, and modify eqn 11.17 for a weighted regression by
%           substituting Sw for n, Swx for Sx, Swy for Sy, Swxy for Sxy, etc. 
%
%     Data are input and output as follows:
%
%         [mw,bw,rw,smw,sbw,xw,yw] = lsqfityw(X,Y,sY)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%             sY    =    estimated uncertainty in y data (vector)
%
%             sy may be measured or calculated:
%                 sY = sqrt(Y), 2% of y, etc.
%             data points are then weighted by:
%                 w = 1 / sY-squared.
%
%             mw    =    slope
%             bw    =    y-intercept
%             rw    =    weighted correlation coefficient
%             smw   =    standard deviation of the slope
%             sbw   =    standard deviation of the y-intercept
%
%     NOTE that the line passes through the weighted centroid: (xw,yw).

function [mw,bw,rw,smw,sbw,xw,yw]=lsqfityw(X,Y,sY)

% Determine the size of the vector

n = length(X);

% Calculate the weighting factors

W = 1 ./ (sY.^2);

% Calculate the sums

Sw = sum(W);
Swx = sum(W .* X);
Swy = sum(W .* Y);
Swx2 = sum(W .* X.^2);
Swxy = sum(W .* X .* Y);
Swy2 = sum(W .* Y.^2);

% Determine the weighted centroid

xw = Swx / Sw;
yw = Swy / Sw;

% Calculate re-used expressions

num = Sw * Swxy - Swx * Swy;
del1 = Sw * Swx2 - Swx^2;
del2 = Sw * Swy2 - Swy^2;

% Calculate mw, bw, rw, smw, and sbw

mw = num / del1;
bw = (Swx2 * Swy - Swx * Swxy) / del1;

rw = num / (sqrt(del1 * del2));

smw = sqrt(Sw / del1);
sbw = sqrt(Swx2 / del1);

return;

function m = meanmiss(x)
%MEANMISS Column matrix means with missing data.
%  MEANMISS(X) returns the means of the columns of X 
%  as a row vector, where missing data in X are encoded as NaN's.
%  Returns the same result as MEAN if there are no NaNs in X.

%m = nan*ones(1,length(x));
%jj = find(sum(finite(x)) ~= 0);
%if isempty(jj)==1, return, end

notvalid = isnan(x);
x(notvalid) = zeros(size(x(notvalid)));
denom = sum(1-notvalid);
if denom==0
	m = nan; return, 
elseif denom==1
	m = x; return
else
	m = sum(x)./denom;
end

return;


% lsqfityz.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
%
% M-file to calculate a "MODEL-1" least squares fit to WEIGHTED x,y-data pairs:
%
%     The line is fit by MINIMIZING the WEIGHTED residuals in Y only.
%
%     The equation of the line is:     Y = mz * X + bz.
%
%     Equations are from Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       for mz and bz, see p. 98, example calculation in Table 6.2;
%       for rz, see p. 199, and modify eqn 11.17 for a weighted regression by
%           substituting Sw for n, Swx for Sx, Swy for Sy, Swxy for Sxy, etc. 
%
%       smz, sbz are adapted from: York (1966) Canad. J. Phys. 44: 1079-1086.
%
%     Data are input and output as follows:
%
%         [mz,bz,rz,smz,sbz,xz,yz] = lsqfityz(X,Y,sY)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%             sY    =    estimated uncertainty in y data (vector)
%
%             sy may be measured or calculated:
%                 sY = sqrt(Y), 2% of y, etc.
%             data points are then weighted by:
%                 w = 1 / sY-squared.
%
%             mz    =    slope
%             bz    =    y-intercept
%             rz    =    weighted correlation coefficient
%             smz   =    standard deviation of the slope
%             sbz   =    standard deviation of the y-intercept
%
%     NOTE that the line passes through the weighted centroid: (xz,yz).

function [mz,bz,rz,smz,sbz,xz,yz]=lsqfityz(X,Y,sY)

% Determine the size of the vector

n = length(X);

% Calculate the weighting factors

W = 1 ./ (sY.^2);

% Calculate the sums

Sw = sum(W);
Swx = sum(W .* X);
Swy = sum(W .* Y);
Swx2 = sum(W .* X.^2);
Swxy = sum(W .* X .* Y);
Swy2 = sum(W .* Y.^2);

% Determine the weighted centroid

xz = Swx / Sw;
yz = Swy / Sw;

% Calculate re-used expressions

num = Sw * Swxy - Swx * Swy;
del1 = Sw * Swx2 - Swx^2;
del2 = Sw * Swy2 - Swy^2;

U = X - xz;
V = Y - yz;
U2 = U.^2;
V2 = V.^2;

% Calculate mw, bw, rw, smw, and sbw

mz = num / del1;
bz = (Swx2 * Swy - Swx * Swxy) / del1;

rz = num / (sqrt(del1 * del2));

sm2 = (1 / (n-2)) * (sum(W .* (((mz * U) - V) .^ 2)) / sum(W .* U2));
smz = sqrt(sm2);
sbz = sqrt(sm2 * (sum(W .* (X.^2)) / Sw));

return;

% lsqfitma.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Sep 05.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The line is fit by MINIMIZING the NORMAL deviates.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the MAJOR AXIS.  All points are given EQUAL
%       weight.  The units and range for X and Y must be the same.
%     Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
%       re-written from Kermack & Haldane (1950) Biometrika 37: 30-41;
%       after a derivation by Pearson (1901) Phil. Mag. V2(6): 559-572.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqfitma(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)
 
function [m,b,r,sm,sb]=lsqfitma(X,Y)

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and other re-used expressions
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;
U = X - xbar;
V = Y - ybar;
 
Suv = sum(U .* V);
Su2 = sum(U .^2);
Sv2 = sum(V .^2);
 
sigx = sqrt(Su2/(n-1));
sigy = sqrt(Sv2/(n-1));
 
% Calculate m, b, r, sm, and sb
 
m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv);
b = ybar - m * xbar;
r = Suv / sqrt(Su2 * Sv2);
 
sm = (m/r) * sqrt((1 - r^2)/n);
sb1 = (sigy - sigx * m)^2;
sb2 = (2 * sigx * sigy) + ((xbar^2 * m * (1 + r))/r^2);
sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);

return;


% lsqfitgm.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The SLOPE of the line is determined by calculating the GEOMETRIC MEAN
%       of the slopes from the regression of Y-on-X and X-on-Y.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the GEOMETRIC MEAN or the REDUCED MAJOR AXIS.
%
%     See Ricker (1973) Linear regressions in Fishery Research, J. Fish.
%       Res. Board Can. 30: 409-434, for the derivation of the geometric
%       mean regression.
%
%     Since no statistical treatment exists for the estimation of the
%       asymmetrical uncertainty limits for the geometric mean slope,
%       I have used the symmetrical limits for a model I regression
%       following Ricker's (1973) treatment.  For ease of computation,
%       equations from Bevington and Robinson (1992) "Data Reduction and
%       Error Analysis for the Physical Sciences, 2nd Ed."  pp: 104, and
%       108-109, were used to calculate the symmetrical limits: sm and sb.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqfitgm(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)
%
%     WARNING:  Both lsqfitx.m and lsqfity.m must be present in a directory
%               that is in your MATLABPATH in order for this algorithm to
%               execute properly.

function [m,b,r,sm,sb]=lsqfitgm(X,Y)

% Determine slope of Y-on-X regression

[my] = lsqfity(X,Y);

% Determine slope of X-on-Y regression

[mx] = lsqfitx(X,Y);

% Calculate geometric mean slope

m = sqrt(my * mx);

if (my < 0) && (mx < 0)
	m = -m;
end

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and means
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;

% Calculate geometric mean intercept

b = ybar - m * xbar;

% Calculate more sums

Sxy = sum(X .* Y);
Sx2 = sum(X.^2);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sx * Sy;
den = n * Sx2 - Sx^2;

% Calculate r, sm, sb and s2

r = sqrt(my / mx);

if (my < 0) && (mx < 0)
	r = -r;
end

diff = Y - b - m .* X;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sb = sqrt(Sx2 * s2 / den);

return;


% lsqbisec.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The SLOPE of the line is determined by calculating the slope of the line
%       that bisects the minor angle between the regression of Y-on-X and X-on-Y.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the LEAST SQUARES BISECTOR.
%
%     See: Sprent and Dolby (1980). The Geometric Mean Functional Relationship.
%       Biometrics 36: 547-550, for the rationale behind this regression.
%
%     Sprent and Dolby (1980) did not present a statistical treatment for the
%       estimation of the uncertainty limits for the least squares bisector
%       slope, or intercept.
%
%     I have used the symmetrical limits for a model I regression following
%       Ricker's (1973) treatment.  For ease of computation, equations from
%       Bevington and Robinson (1992) "Data Reduction and Error Analysis for
%       the Physical Sciences, 2nd Ed."  pp: 104, and 108-109, were used to
%       calculate the symmetrical limits: sm and sb.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqbisec(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)
%
%     WARNING:  Both lsqfitx.m and lsqfity.m must be present in a directory
%               that is in your MATLABPATH in order for this algorithm to
%               execute properly.
 
function [m,b,r,sm,sb]=lsqbisec(X,Y)

% Determine slope of Y-on-X regression

[my] = lsqfity(X,Y);

% Determine slope of X-on-Y regression

[mx] = lsqfitx(X,Y);

% Calculate the least squares bisector slope

theta = (atan(my) + atan(mx)) / 2;
m = tan(theta);

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and means
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;

% Calculate the least squares bisector intercept

b = ybar - m * xbar;

% Calculate more sums

Sxy = sum(X .* Y);
Sx2 = sum(X.^2);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sx * Sy;
den = n * Sx2 - Sx^2;

% Calculate r, sm, sb and s2

r = sqrt(my / mx);

if (my < 0) && (mx < 0)
	r = -r;
end

diff = Y - b - m .* X;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sb = sqrt(Sx2 * s2 / den);

return;

% lsqfitx.m                                      by:  Edward T Peltzer, MBARI
%                                                revised:  2007 Apr 28.
%
% M-file to calculate a "MODEL-1" least squares fit.
%
%     The line is fit by MINIMIZING the residuals in X only.
%
%     The equation of the line is:     Y = mx * X + bx.
%
%     Equations are modified from those in Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       pp: 104, 108-109, 199.
%
%     Data are input and output as follows:
%
%         [mx,bx,rx,smx,sbx] = lsqfitx(X,Y)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%
%             mx     =    slope
%             bx     =    y-intercept
%             rx     =    correlation coefficient
%             smx    =    standard deviation of the slope
%             sbx    =    standard deviation of the y-intercept

function [mx,bx,rx,smx,sbx]=lsqfitx(X,Y)

% Determine the size of the vector

n = length(X);

% Calculate the sums

Sx = sum(X);
Sy = sum(Y);
Sx2 = sum(X.^2);
Sxy = sum(X.*Y);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sy * Sx;
den = n * Sy2 - Sy^2;

% Calculate m, a, rx, s2, sm, and sb

mxi = num / den;
a = (Sy2 * Sx - Sy * Sxy) / den;
rx = num / (sqrt(den) * sqrt(n * Sx2 - Sx^2));

diff = X - a - mxi .* Y;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sa = sqrt(Sy2 * s2 / den);

% Transpose coefficients

mx = 1 / mxi;
bx = -a / mxi;

smx = mx * sm / mxi;
sbx = abs(sa / mxi);

return;



