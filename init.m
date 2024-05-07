close; format long
for i=1:7
input=load([dir{3},'init/input.data.', Ray{i}]);
% outInput=input(:,1:5);



bet1=3/40.0;
omega=input(:,5)./input(:,4)/0.09;

omega(1)=60*nu(1,i)/(bet1*(input(2,1)-input(1,1))^2);
omega(end)=60*nu(1,i)/(bet1*(input(N(i),1)-input(N(i)-1,1))^2);

%omega(1)=100.;
%omega(N)=100.;

outInput=[input(:,1:4), omega];

% input(:,1:4);


save(['init/input.data.',Ray{i}],'outInput','-ASCII')
% close all
% N=201;
% % Ra = 5.4e+5;
% % g = 9.81 ; % gravity acceleration 
% % N=251;
% % Ra = 1.0e+9;
% Pr = 0.709; % prandtl number
% Gr=Ra/Pr;
% % beta=1.0/Gr/g;
% nu=1.0/Gr;
end