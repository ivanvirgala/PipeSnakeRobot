function res= objective(t,param,optParametre,x0)
    t=0:param.dt:20;
    [T,X] = ode45(@(t,y)dynamicModelOptimization(t,y,param,optParametre),t,x0);
    res = abs(X(length(X(:,param.N+1)),param.N+1) - X(2,param.N+1)) 
    optParametre
    disp('-------------')    
end