function []=font(sz)

%font(sz)
%
%set font size of a figure

if nargin==0,
sz=9;
end;

set(gca,'fontsize',sz);
