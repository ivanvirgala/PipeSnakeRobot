cas=0;
dt = 0.01;
for i=1:length(X)
    distanceV(i,1) = sqrt(X(i,16)^2 + X(i,17)^2);
    distanceV(i,2) = cas;
    cas = cas + dt;
end
hold on
plot(distanceV(:,2), distanceV(:,1))