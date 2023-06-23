function res= objective(t,param,optParametre,x0)
    t=0:param.dt:60;
    [T,X] = ode45(@(t,y)dynamicModelOptimization(t,y,param,optParametre),t,x0);  % -  kvoli hladaniu maxima, inak by sa hladalo minimum
    res = abs(X(length(X(:,param.N+1)),param.N+1) - X(2,param.N+1)) 
    optParametre
    disp('-------------')
    %assignin('base', 'x_saved', res);
    %res = abs(X(length(X(:,param.N+1)),param.N+1) - X(2,param.N+1))
    
end