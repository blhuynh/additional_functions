function gzoom(n, f)
% gzoom.m
% Michael R. Gustafson II (mrg@duke.edu)
% Uploaded to pundit.pratt.duke.edu on September 6, 2008
%
% Function to determine the size of a figure window and increase the size
% GZOOM will increase figure 1 by a factor of 
%       0.1 of the original size in each direction
% GZOOM(r) will increase figure 1 by a factor of
%       r of the original size in each direction   
% GZOOM(r, f) will increase figure f by a factor of 
%       r of the original size in each direction
 
if nargin==0
    f=gca; n=0.1;
elseif nargin==1
    f=gca;
end
v=axis;
xr = (v(2)-v(1));
yr = (v(4)-v(3));
axis(f, [v + n*[-xr xr -yr yr]]);