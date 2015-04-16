function [v,f,xgrid,ygrid]= probability_contours(x,fracs)
%[v,f,xgrid,ygrid]=prob_contour(x,fracs);  x is n by 2 matrix of data,
%fracs is optional fractions where you (approximately) want the contours,
%for example [0.2, 0.4, 0.6, 0.8] for four contours at the 20th, 40th, 60th
%and 80th percentiles.  without this second argument, the default is
%0.1:0.1:0.9.  v is a vector listing the density levels selected for the
%contours; f is the matrix of density estimates, and xgrid and ygrid are
%the meshgrid for the x- and y- coordinates. to make a contour plot, you
%would then use contour(xgrid,ygrid,f,v)

if nargin<2
fracs=0.1:0.1:0.9;  %the percentages at which you want to plot contourss
end

[f,xgrid,ygrid]=densmat(x,50);
f=f';

dA=prod(([xgrid(end) ygrid(end)]-[xgrid(1) ygrid(1)])/127);

fdA=f(:)*dA;
maxf=max(fdA(:));
fints=linspace(0,maxf,1000);
s=0;
v=zeros(size(fracs));
probs=zeros(size(fracs));
j=1;
for i=1:length(fracs)    
    while s<fracs(i) && j <= length(fints)
    fi=fints(j);
    s=sum(fdA(fdA<fi));
    j=j+1;
    end
   v(i)=fi/dA;  %these will be the density levels corresponding to the probabilities listed in fracs
   probs(i)=s; 
end

%you can check the accuracy by comparing probs to fracs. adding more points
%to fints or changing its segmentation will help

%now, plot the contours at specified density levels v
% contour(xgrid',ygrid',f',v,'k')

function [fmat,xm,ym]=densmat(x,perc,M)

if nargin<3
    M=128;
end

if nargin<2
    perc=0;
end

  %number of gridpoints in each dim

%GridAssign=bin2(x(:,1),x(:,2),M);

n=size(x,1);  %number of data points
d=size(x,2);  %number of dimensions 
MM = M^d; %number of total gridpoints

mins=min(x,[],1);
maxs=max(x,[],1);
Diff=maxs-mins;
mins=mins-perc/100*Diff;
maxs=maxs+perc/100*Diff;

Delta = 1/(M-1)*(maxs-mins);  %vector of distances between neighbor grid points in each dimension


ye=zeros(d,M);
multby=zeros(1,d);  %used in coord transfrom from m to k
pointLL=zeros(n,d);  %this will be the "lower left" gridpoint to each data point
for i = 1:d
    ye(i,:) = linspace(mins(i),maxs(i),M);
    multby(i)=M^(i-1);
    pointLL(:,i)=floor((x(:,i)-mins(i))./Delta(i)) + 1;
end
pointLL(pointLL==M)=M-1;  %this avoids going over grid boundary

%% assign each data point to its closest grid point
[xgrid,ygrid]=meshgrid(ye(1,:),ye(2,:));
z=reshape(1:MM,M,M);
GridAssign=interp2(xgrid,ygrid,z',x(:,1),x(:,2),'nearest');  %this associates each data point with its nearest grid point

%% compute w
Deltmat=repmat(Delta,n,1);
shape=M*ones(1,d);

wmat=zeros(M,M);
for i=0:1  %number of neighboring gridpoints in d dimensions
    for j=0:1
        pointm=pointLL+repmat([j i],n,1);  %indices of ith neighboring gridpoints
        pointy=zeros(n,d);
        for k=1:d
            pointy(:,k)=ye(k,pointm(:,k));  %y-values of ith neighboring gridpoints
        end
        W=prod(1-(abs(x-pointy)./Deltmat),2);  %contribution to w from ith neighboring gridpoint from each datapoint
        wmat=wmat+accumarray(pointm,W,shape);  %sums contributions for ith gridpoint over data points and adds to wmat
    end
end

%% compute f, sig, df and A
n6 = n^(-1/6);
h=zeros(1,d);
Z=zeros(1,d);
Zin=cell(1,d);
for i =1:d
    h(i) = std(x(:,i))*n6;
    Z(i) = min(floor(4*h(i)/Delta(i)), M-1);
    Zin{i}=-Z(i):Z(i);
end

phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);

[L{1},L{2}]=meshgrid(Zin{1},Zin{2});

Phix=phi(L{1}*Delta(1)./h(1))./h(1);
Phiy=phi(L{2}*Delta(2)./h(2))./h(2);

Phimat = (Phix.*Phiy)';   %matrix of Phi for inputting into convn

fmat = 1/n*conv2(wmat,Phimat,'same');  %d-dim matrix of estimated densities
[xm,ym]=meshgrid(ye(1,:),ye(2,:));

% if nargout==0
% contour(xm,ym,fmat',20,'k')
% end