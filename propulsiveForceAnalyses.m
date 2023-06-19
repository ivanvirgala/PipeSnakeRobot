clear
clc


N = 15;
% Viscous
k = ones(1,N)';
ct = 0.015;
cn = 0.03;
ut = 0.15;
un = 0.30;
m = 0.406;
g = 9.81;
xc = 1;
yc = 1;

i = 1;
for q1=-90:1:90
    fx(i) = ct*cosd(q1)*cosd(q1) + cn*sind(q1)*sind(q1);
    fy(i) = (ct-cn)*sind(q1)*cosd(q1);
    i = i +1;
end
%plot(fx)
%hold on
plot(fy)
%legend('fx','fy')
grid on
%{
% Coulomb
i = 1;
for q1=-90:1:90
    fx(i) = sind(q1)*g*m*un*sign(cosd(q1)*yc - sind(q1)*xc) - cosd(q1)*g*m*ut*sign(cosd(q1)*xc + sind(q1)*yc);
    fy(i) = - cosd(q1)*g*m*un*sign(cosd(q1)*yc - sind(q1)*xc) - sind(q1)*g*m*ut*sign(cosd(q1)*xc + sind(q1)*yc);
    i = i +1;
end
plot(fx)
hold on
plot(fy)
legend('fx','fy')
%}