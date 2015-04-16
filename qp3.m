function qp3(data)

if nargin > 1; error('One input only'); end
if ~any(size(data)==3); error('Must have 3 columns or rows'); end
if size(data,1)==3; data=data'; end

%% Plot Stuff
figure; hold on
plot3(data(:,1),data(:,2),data(:,3))
view(3)
grid on

xlabel('x')
ylabel('y')
zlabel('z')
setfont

end

