clear
clc

N = 15;
k = ones(1,N)';
ct = 0.015;
cn = 0.03;


i = 1;
for q1=-90:1:90
    fx(i) = ct*cosd(q1)*cosd(q1) + cn*sind(q1)*sind(q1);
    fy(i) = (ct-cn)*sind(q1)*cosd(q1);
    i = i +1;
end
plot(fx)
hold on
plot(fy)
legend('fx','fy')
grid on
