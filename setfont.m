function setfont(varargin)

if nargin==0
    fontSize=14;
elseif nargin==1
    fontSize=varargin{1};
end

set(gca,'FontSize',fontSize)
figureHandle=gcf;
set(findall(figureHandle,'type','text'),'fontSize',fontSize)