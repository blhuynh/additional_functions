function [yvalid, nvalid, phandle]=ConvPro(x, h, nxi, nhi, xtrunc, htrunc, pflag)
% CONVPRO   Dr. G's Convolver Pro v. 1.1.2
% Perform discrete convolution on finite duration and/or 
% right-sided signals with truncated representations.
%
% [YVALID, NVALID] = CONVPRO(X, H) convolves vectors X and H and returns
% the result as well as the time base for the result.
% 
% [YVALID, NVALID] = CONVPRO(X, H, NXI, NHI) is the same as above, only
% the first entry in X is assumed to be at time n=NXI and
% the first entry in H is assumed to be at time n=NHI
%
% [YVALID, NVALID] = CONVPRO(X, H, NXI, NHI, XTRUNC, HTRUNC) is the same
% as above, only now the program will check the values of XTRUNC and
% HTRUNC. If one of the truncation variables is 1, the program will
% assume the relevant vector is a finite model for a signal of 
% infinite duration and truncate the output vector and timebase to
% contain only those which are valid given the length of the model.
%
% [YVALID, NVALID, PHANDLE] = CONVPRO(X, H, NXI, NHI, XTRUNC, HTRUNC, PFLAG) is the
% same as the above, only now the program will check the value of
% PFLAG.  If it is non-zero, a plot of the input, impulse response, and
% output will also be generated and handles to the three subplots will
% be returned to PHANDLE

% Last updated 9/27/2013 by Michael R. Gustafson II
% Copyright 2007-2013
% Revision: 1.0.0 - 2/10/2007 - initial rollout
% Revision: 1.0.1 - 2/14/2010 - better comments
% Revision: 1.1.1 - 2/15/2010 - passes print handles
% Revision: 1.1.2 - 9/27/2013 - updated variable names
%
% To Do: left-sided signal flags

%% Error checking
if nargin<2, error('Must have 2 signals!'); end
if nargin<7 & nargout==3, error('Cannot return handles if not plotting!'); end

%% Default cases
% Default starting times are 0
if nargin==2, nxi = 0; end
if nargin<=3, nhi = 0; end
% Default truncation flags are off
if nargin<=4, xtrunc = 0; end
if nargin<=5, htrunc = 0; end
% Default plotting flag is off
if nargin<=6, pflag=0; end

%% Convolve x and h
y = conv(x, h);

%% Time-base validation and limitation
% Determine time bases of signals
nx = nxi + (0:(length(x)-1));
nh = nhi + (0:(length(h)-1));
ny = (nxi+nhi)+(0:length([x h 1]));
% Determine valid range of output
xvalid = length(x)+realmax*(xtrunc==0);
hvalid = length(h)+realmax*(htrunc==0);
NT = min(min(xvalid, hvalid), length(y));
yvalid = y(1:NT);
nvalid = ny(1:NT);

%% Print if requested
if pflag==1
    % Determine plot limits
    nleft  = min([nx nh nvalid])-1;  nright = max([nx nh nvalid])+1;
    % Plot three signals
    phandle(1)=subplot(3,1,1); stem(nx, x, 'k');
    axis([nleft nright get(gca, 'YLim')+[-.1 .1]*max(abs(x))])
    xlabel('n'); ylabel('x[n]'); title('Input Signal: x[n]');
    phandle(2)=subplot(3,1,2);
    stem(nh, h, 'k');
    axis([nleft nright get(gca, 'YLim')+[-.1 .1]*max(abs(h))])
    xlabel('n'); ylabel('h[n]'); title('Impulse Response: h[n]');
    phandle(3)=subplot(3,1,3);
    stem(nvalid, yvalid, 'k')
    hold on
    stem(nvalid(end)+1:nright-1, 0*(nvalid(end)+1:nright-1), 'rx')
    stem(nleft+1:nvalid(1)-1, 0*(nleft+1:nvalid(1)-1), 'bx')
    axis([nleft nright get(gca, 'YLim')+[-.1 .1]*max(abs(yvalid))])
    xlabel('n'); ylabel('y[n]');
    title(sprintf('Output Signal: y[n] valid for n=[%0.0f, %0.0f]',...
        nvalid(1), nvalid(end)));
    hold off
end