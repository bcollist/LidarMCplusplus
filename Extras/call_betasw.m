%% Calculate the VSF of Seawater

lambda = [400:1:700];
Tc = 20;
theta = [0:90/454:180]'
S = 35;

[betasw,beta90sw,bsw]= betasw_ZHH2009(lambda,Tc,theta,S);

m11 = betasw;

m12n = -sind(theta).^2 ./ (1+cosd(theta).^2)
m33n = 2.*cosd(theta)./(cosd(theta).^2 +1)

m12 = m12n .* m11;
m33 = m33n .* m11;

mExport = [m11 m12 m33];
%% Export File
csvwrite('/Users/Brian/Documents/C++/LidarMCplusplus/seawaterVSFZHHM11.csv',m11)
csvwrite('/Users/Brian/Documents/C++/LidarMCplusplus/seawaterVSFZHHM12.csv',m12)
csvwrite('/Users/Brian/Documents/C++/LidarMCplusplus/seawaterVSFZHHM33.csv',m33)









