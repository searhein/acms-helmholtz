function a=set_diffusion_coeffcient(p,t,n)
cv1= 1 ; % first value
cv2= 12; % second value within perturbation
a=zeros(1,size(t,2))+ cv1;
p1=p(:,t(1,:));
p2=p(:,t(2,:));
p3=p(:,t(3,:));
m=(p1+p2+p3)/3; % midpoints
x=m(1,:); % in [0,1];
y=m(2,:); % in [0,1];
X=n*x;  % in [0,n];
Y=n*y;  % in [0,n];
i=floor(X); % integer part of X
j=floor(Y); % integer part of Y
xx=X-i; % noninteger part of X in [0,1]
yy=Y-j; % noninteger part of Y in [0,1]
ind=(xx>1/4) & (xx<3/4) & (yy>1/4) & (yy<3/4); % for rectangular inclusion
%ind=((xx-1/2).^2+(yy-1/2).^2 <1/16); % for circular inclusion
a(ind)=cv2;
% a=exp(-10*((xx-1/2).^2+(yy-1/2).^2)); % smooth diffusion coefficient
end