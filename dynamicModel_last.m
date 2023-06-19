function [xDot] = dynamicModel(t,x,param)
    %% Basic parameters
    ut      = param.ut;
    un      = param.un;
    ct      = param.ct;
    cn      = param.cn;
    N       = param.N;
    kp      = param.kp;
    kd      = param.kd;
    l       = param.l;
    d       = param.d;   
    priemer = param.priemer;
    priemerI= param.priemerInfluence;
    m       = param.m;
    g       = param.g;
    alfa    = param.alfa;
    omega   = param.omega;
    delta   = param.delta;
    offset  = param.offset;
    tlmic   = param.tlmic;
    pruzina = param.pruzina;
    viskozne= param.trenie;
    kontakt = param.kontakt;
    utPipe  = param.utPipe;
    ctPipe  = param.ctPipe;
    Erub    = param.Erub;
    vrub    = param.vrub;
    umax    = param.umax;
    qmax    = param.qmax;
    minLinkVel=param.minLinkVel;

    d = priemer - 2*l;
    I = (m*(2*l)^2)/3;  
    order=1;
    %% Auxiliary matrices
    for i=1:N-1
        for ii = 1:(N)
            if(i==ii)
                Alil(i,ii) = 1;
            end
            if(ii==i+1)
                Alil(i,ii) = 1;
            end
        end
    end
    
    for i=1:N-1
        for ii = 1:(N)
            if(i==ii)
                Dlil(i,ii) = 1;
            end
            if(ii==i+1)
                Dlil(i,ii) = -1;
            end
        end
    end
    
    Vlil = Alil'*inv(Dlil*Dlil')*Alil;
    Klil = Alil'*inv(Dlil*Dlil')*Dlil;
        
    for i=1:N
        for ii = 1:(N+1)
            if(i==ii)
                J(i,ii) = -1;
            end
            if(ii==i+1)
                J(i,ii) = 1;
            end
        end
    end
    K = zeros(N,N+1);
    for i=1:N
        for ii = 1:N+1
            if(i==ii)
                K(i,ii) = 1;
            end
        end
    end
    
    for i=1:N
        for ii = 1:(N-1)
            if(i==ii)
                J2(i,ii) = 1;
            end
            if(i==ii+1)
                J2(i,ii) = 1;
            end
        end
    end
    
    for i=1:N
        for ii = 1:(N-1)
            if(i==ii)
                J4(i,ii) = 1;
            end
            if(i==ii+1)
                J4(i,ii) = -1;
            end
        end
    end
    
    HH = ones(N,N);
    HH = triu(HH);
        
    J3 = -J2;
    J1 = -J4;
    NN = J*pinv(K);
    T = abs(NN);
    
    e = ones(1,N+1)'; e(length(e)) = 0;
    k = ones(1,N)';   
    j = zeros(1,N)';
    j(N) = -1;   

    %% Kinematika
    % state vector
    fi    = x(1:N);
    p     = x(N+1:N+2);
    fiDot = x(N+3:2*N+2);
    pDot  = x(2*N+3:2*N+4);

    theta = HH*fi;
    thetaDot = HH*fiDot;

    for oo=1:N
        sss(oo,1) = sin(theta(oo));
        ccc(oo,1) = cos(theta(oo));
        sgn(oo,1) = sign(theta(oo));
        dThetaqSqared(oo,1) = thetaDot(oo)^2;
    end
    
    for oo=1:N
        for ooo=1:N
            if(oo==ooo)
                Cm(oo,ooo) = cos(theta(ooo));
                Sm(oo,ooo) = sin(theta(ooo));
            end
        end
    end

    A = J*J';
    B = (1/N)*j;
    C = (1/N)*j';
    D = (1/N);

    DD = inv(D-C*inv(A)*B);
    AA = inv(A)+inv(A)*B*DD*C*inv(A);
    BB = -inv(A)*B*DD;
    CC = -DD*C*inv(A);

    Hinv = [J' (1/N)*e]*[AA BB;CC DD];
    
    X = Hinv*[2*l*ccc;p(1) - (l/N)*k'*ccc];
    Y = Hinv*[2*l*sss;p(2) - (l/N)*k'*sss];
    
    dX = -J'*AA*2*l*Sm*thetaDot - (1/N)*e*CC*2*l*Sm*thetaDot + J'*BB*pDot(1) + J'*BB*(l/N)*k'*Sm*thetaDot + (1/N)*e*DD*pDot(1) + (1/N)*e*DD*(l/N)*k'*Sm*thetaDot;
    dY =  J'*AA*2*l*Cm*thetaDot + (1/N)*e*CC*2*l*Cm*thetaDot + J'*BB*pDot(2) - J'*BB*(l/N)*k'*Cm*thetaDot + (1/N)*e*DD*pDot(2) - (1/N)*e*DD*(l/N)*k'*Cm*thetaDot;
    
    Xc = K*X + l*ccc;
    Yc = K*Y + l*sss;
    
    dXc = K*dX - l*Sm*thetaDot;
    dYc = K*dY + l*Cm*thetaDot;
    
    %% Friction

    % Viscous
    if(viskozne==1)
        fr = -[ct*(Cm.^2)+cn*(Sm.^2), (ct-cn)*Sm*Cm;(ct-cn)*Sm*Cm, ct*(Sm.^2)+cn*(Cm.^2)]*[dXc;dYc];
    else    
    % Coulomb
        fr = -m*g*[ut*Cm, -un*Sm;ut*Sm, un*Cm]*sign([Cm, Sm;-Sm, Cm]*[dXc;dYc]);
    end
    
    %% Contact
    % ============================== Kontakt ==============================
    if kontakt == 1
        for mm=1:N
            fcn(mm,1) = 0;
            fctBool(mm,1) = 0;
            fct(mm,1) = 0;
        end
        
        for mm=1:N  

            if((Yc(mm)>=(d/2)))
                fcn(mm,1) = -(sqrt((16*Erub*2*l)/(9*(1-vrub^2))))*((abs(Yc(mm))-(d/2))^(3/2));
                if(abs(dXc(mm))>minLinkVel)
                    fctBool(mm,1) = 1;
                    if(viskozne == 1)
                        fct(mm,1) = fctBool(mm)*ctPipe*sign(-dXc(mm));
                    else
                        fct(mm,1) = fctBool(mm)*(abs(fcn(mm))*utPipe*sign(-dXc(mm)));
                    end
                else
                    fctBool(mm,1) = 0;
                    fct(mm,1) = 0;
                end
            elseif((Yc(mm)<=(-d/2)))
                fcn(mm,1) = (sqrt((16*Erub*2*l)/(9*(1-vrub^2))))*((abs(Yc(mm))-(d/2))^(3/2));
                if(abs(dXc(mm))>minLinkVel)
                    fctBool(mm,1) = 1;
                    if(viskozne == 1)
                        fct(mm,1) = fctBool(mm)*ctPipe*sign(-dXc(mm));
                    else
                        fct(mm,1) = fctBool(mm)*(abs(fcn(mm))*utPipe*sign(-dXc(mm)));
                    end
                else
                    fctBool(mm,1) = 0;
                    fct(mm,1) = 0;
                end
            else
                fcn(mm,1) = 0;
                fctBool(mm,1) = 0;
                fct(mm,1) = 0;
            end
            
        end

        fcontact = [fct;fcn];
    else
        fcontact = zeros(N*2,1);
        fct = zeros(N,1);
        fcn = zeros(N,1);
        fctBool = zeros(N,1);
    end

    ground = fcontact + fr;

    %% Model

    for i=1:N-1
        fi_required         = alfa*sin((omega*t+(i-1)*delta)) + offset;
        fiDot_required      = alfa*omega*cos((omega*t+(i-1)*delta));
        fiDotDot_required   = -alfa*omega*omega*sin((omega*t+(i-1)*delta));
        %u(i,1)              = kp*(fi_required-fi(i))-kd*fiDot(i);                  %control law
        u(i,1)              = fiDotDot_required + kp*(fi_required - fi(i)) + kd*(fiDot_required - fiDot(i));
        if(u(i,1)>umax)
            u(i,1) = umax;
        elseif(u(i,1)<-umax)
            u(i,1) = -umax;
        end
    end

    M = I*eye(N) + m*l*l*Sm*Vlil*Sm + m*l*l*Cm*Vlil*Cm;
    W = m*l*l*Sm*Vlil*Cm - m*l*l*Cm*Vlil*Sm;

    fi2 = fi;
    fi2(N,1) = theta(N);
    fi2Dot = fiDot;
    fi2Dot(N,1) = thetaDot(N);

    MMM = [HH'*M*HH zeros(N,2);zeros(2,N) N*m*eye(2,2)];
    WWW = [HH'*W*diag(HH*fiDot)*HH*fiDot;zeros(2,1)];
    GGG = [-l*HH'*Sm*Klil l*HH'*Cm*Klil;-k' zeros(1,N);zeros(1,N) -k'];
    BBB = [eye(N-1,N-1);zeros(3,N-1)];

    M11 = MMM(1:N-1,1:N-1);
    M12 = MMM(1:N-1,N:N+2);
    M21 = MMM(N:N+2,1:N-1);
    M22 = MMM(N:N+2,N:N+2);
    W1  = WWW(1:N-1);
    W2  = WWW(N:N+2);
    G1  = GGG(1:N-1,1:2*N);
    G2  = GGG(N:N+2,1:2*N);
    Aq  = -inv(M22)*(W2 + G2*fr + G2*fcontact);
    Bq  = -inv(M22)*M21;
    xDot = [fiDot;pDot;u;Aq+Bq*u;fcontact];
    %xDot = [fiDot;pDot;u;Aq+Bq*u];
    
    %{
    fiKryt = 0:0.1:2*pi;
    xKryt = (param.l)*cos(fiKryt);
    yKryt = (param.l)*sin(fiKryt);

    for k=1:N
        %h(k) = line([xLink(k,i) xLink(k,i)],[yLink(k,i) yLink(k,i)],'color','red','LineWidth',2);
        hold on;
        %h(k) = line([result_X(k,i) (result_X(k,i)+l*cosd(result_theta(k,i)))],[result_Y(k,i) (result_Y(k,i)+l*sind(result_theta(k,i)))],'color','black');
        %hold on
        if(Yc(k)==0) %fcn(k)
            cmFct(k) = line([Xc(k) (Xc(k)-0.1)],[Yc(k) (Yc(k))],'color','white','LineWidth',3);%X(k,(2*param.N+4)+k))
            hold on
            cmFcn(k) = line([Xc(k) (Xc(k))],[Yc(k) (Yc(k)+0.1)],'color','white','LineWidth',3);%+X(k,(2*param.N+4)+k))
            hold on
        else


            cmFct(k) = line([Xc(k) (Xc(k)-0.1)],[Yc(k) (Yc(k))],'color','green','LineWidth',3);%X(k,(2*param.N+4)+k))
            hold on
            cmFcn(k) = line([Xc(k) (Xc(k))],[Yc(k) (Yc(k)+0.1)],'color','green','LineWidth',3);%+X(k,(2*param.N+4)+k))
            hold on
        end
        if(k==1)
            kryt(k) = plot(Xc(k)+xKryt,Yc(k)+yKryt,'color','blue','LineWidth',2);
        else
            kryt(k) = plot(Xc(k)+xKryt,Yc(k)+yKryt,'color','blue','LineWidth',2);
        end
        hold on;
    end
    grid on;
    
    pause(.001);
    
    if(i<length(X(:,1)))
        for lm=1:N
            %delete(p(lm));
            %delete(h(lm));
            %delete(m(lm));
            delete(cmFct(lm));delete(cmFcn(lm));
            %delete(silX);delete(silY);
            delete(kryt(lm));
        end
    end
    %}
end
    


