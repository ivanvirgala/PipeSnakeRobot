%% Animation of snake robot motion
figure('units','normalized','outerposition',[0 0 1 1]);
fiKryt = 0:0.1:2*pi;
xKryt = (param.l)*cos(fiKryt);
yKryt = (param.l)*sin(fiKryt);
if((param.kontakt == 1) & (param.dimensionPlot3D == 0))
    % Inside of the pipe
    line([-dlzkaPotrubia dlzkaPotrubia],[0.025+param.d/2 0.025+param.d/2],'color','black');
    hold on
    line([-dlzkaPotrubia dlzkaPotrubia],[-0.025-param.d/2 -0.025-param.d/2],'color','black');
    hold on
    
    % Outside of the pipe
    line([-dlzkaPotrubia dlzkaPotrubia],[0.025+param.priemer/2 0.025+param.priemer/2],'color','black', 'linewidth', 5); %10je iba kvoli zozbazeniu aby segment robota neprekracoval potrubie
    hold on
    line([-dlzkaPotrubia dlzkaPotrubia],[-0.025-param.priemer/2 -0.025-param.priemer/2],'color','black', 'linewidth', 5);
    hold on
elseif((param.kontakt == 1) & (param.dimensionPlot3D == 1))
    % 3D pipe
    plotcube([3,0.1,param.d],[-2,-param.priemer-0.025,0],0.3,[0 1 0]);
    hold on
    plotcube([3,0.1,param.d],[-2,param.priemer+0.025,0],0.3,[0 1 0]);
    hold on
    plotcube([3,2*param.priemer-0.05,0.01],[-2,-param.priemer+0.075,0],0.3,[0 1 0]);
    hold on
    camlight;
    lighting phong;
    hold on
end

% Outside of the pipe
%{
line([-dlzkaPotrubia dlzkaPotrubia],[param.priemerInfluence/2 param.priemerInfluence/2],'color','green')
hold on
line([-dlzkaPotrubia dlzkaPotrubia],[-param.priemerInfluence/2 -param.priemerInfluence/2],'color','green')
hold on
%}

% Auxiliary matrices
HH = ones(param.N,param.N);
HH = triu(HH);

for i=1:param.N
    for ii = 1:(param.N+1)
        if(i==ii)
            J(i,ii) = -1;
        end
        if(ii==i+1)
            J(i,ii) = 1;
        end
    end
end

K = zeros(param.N,param.N+1);
for i=1:param.N
    for ii = 1:param.N+1
        if(i==ii)
            K(i,ii) = 1;
        end
    end
end

k = ones(1,param.N)';   
j = zeros(1,param.N)';
j(param.N) = -1;
e = ones(1,param.N+1)'; e(length(e)) = 0;

for i=1:length(X(:,1))
    thetaa(:,i) = HH*X(i,1:param.N)';
end

for oo=1:length(thetaa(1,:))
    sss(:,oo) = sin(thetaa(:,oo));
    ccc(:,oo) = cos(thetaa(:,oo));
end



A = J*J';
B = (1/param.N)*j;
C = (1/param.N)*j';
D = (1/param.N);

DD = inv(D-C*inv(A)*B);
AA = inv(A)+inv(A)*B*DD*C*inv(A);
BB = -inv(A)*B*DD;
CC = -DD*C*inv(A);

Hinv = [J' (1/param.N)*e]*[AA BB;CC DD];

for i=1:length(ccc(1,:))
    xLink(:,i) = Hinv*[2*param.l*ccc(:,i);X(i,param.N+1) - (param.l/param.N)*k'*ccc(:,i)];
    yLink(:,i) = Hinv*[2*param.l*sss(:,i);X(i,param.N+2) - (param.l/param.N)*k'*sss(:,i)];
    xcLink(:,i) = K*xLink(:,i) + param.l*ccc(:,i);
    ycLink(:,i) = K*yLink(:,i) + param.l*sss(:,i);
end

saveFig = 0;
fig=figure(1);
row = 3;
column = 3;
order = 1;
del = 0;

% sphere
[XXX,YYY,ZZZ] = sphere;
xx = XXX * param.d/2;
yy = YYY * param.d/2;
zz = ZZZ * param.d/2;

for i=1:pocetSnimkov:length(X(:,1)) 
    %subplot(row,column,order) ............................................odstranit
    if(param.dimensionPlot3D == 1)
        
        % ============================== 3D ===============================
        for k=1:param.N
            cm(k) = surf(xx+xcLink(k,i),yy+ycLink(k,i),zz+param.d/2,'FaceColor','blue','EdgeColor','none');
            hold on
        end

        saveFig = saveFig + 1;
        modulo = mod(saveFig,7);
saveas(fig,sprintf('Fig/3D_FIG_%d.jpg',saveFig));

        if(param.resultsShow == 1)
            if(modulo == 0)
                del = 1;
                
                hold on
                %{
                for j=1:100:i
                    cm(k) = surf(xx+xcLink(k,i),yy+ycLink(k,i),zz+param.d/2,'FaceColor','red','EdgeColor','none');
                    hold on
                end
                %}
                set(gca,'fontsize',14);
                %xlabel('X-axis (m)')
                %ylabel('Y-axis (m)')
                %str = sprintf("Time: %0.2f (s)",(i/100));
                %title(str)
                saveas(fig,sprintf('Fig/3D_FIG_%d.jpg',saveFig));
                order = order + 1;
                delete(cm(1));
            end 
        end

        % ---------------------------- 3D pipe ----------------------------
        plotcube([3,0.1,param.d],[-2,-param.priemer-0.025,0],0.3,[0 1 0]);
        hold on
        plotcube([3,0.1,param.d],[-2,param.priemer+0.025,0],0.3,[0 1 0]);
        hold on
        plotcube([3,2*param.priemer-0.05,0.01],[-2,-param.priemer+0.075,0],0.3,[0 1 0]);
        hold on
        camlight;
        lighting phong;
        hold on
        %------------------------------------------------------------------
        axis equal
        %axis([-2 2 -0.5 0.5 0 0.5])
        pause(0.1);


        
        %if(i<length(X(:,1)))
        if((i<length(X(:,1)))& (del == 0))% & (del == 0)
            for lm=1:param.N
                delete(cm(lm));
            end
        end

        del = 0;

    else
        % ============================== 2D ===============================
        cm = plot(X(i,param.N+1),X(i,param.N+2),'.','color','red');
        eval(['xlabel(',(num2str(round(T(i,1)))),')']);
        for k=1:param.N
            h(k) = line([xLink(k,i) xLink(k,i)],[yLink(k,i) yLink(k,i)],'color','red','LineWidth',2);
            hold on
            
            if(((ycLink(k,i)>(param.d/2)) | (ycLink(k,i)<(-param.d/2))) & (param.kontakt == 1)) 
                %kryt(k) =
                %plot(xcLink(k,i)+xKryt,ycLink(k,i)+yKryt,'color','green','LineWidth',2);
                %pre kontakt
                kryt(k) = plot(xcLink(k,i)+xKryt,ycLink(k,i)+yKryt,'color','blue','LineWidth',2);
            else
                kryt(k) = plot(xcLink(k,i)+xKryt,ycLink(k,i)+yKryt,'color','blue','LineWidth',2);
            end
              
            hold on;
        end


        saveFig = saveFig + 1;
        modulo = mod(saveFig,6);
        if param.kontakt == 1
             % Inside of the pipe
            line([-dlzkaPotrubia dlzkaPotrubia],[0.03+param.d/2 0.03+param.d/2],'color','black');
            hold on
            line([-dlzkaPotrubia dlzkaPotrubia],[-0.03-param.d/2 -0.03-param.d/2],'color','black');
            hold on
            
            % Outside of the pipe
            line([-dlzkaPotrubia dlzkaPotrubia],[0.03+param.priemer/2 0.03+param.priemer/2],'color','black', 'linewidth', 3); %10je iba kvoli zozbazeniu aby segment robota neprekracoval potrubie
            hold on
            line([-dlzkaPotrubia dlzkaPotrubia],[-0.03-param.priemer/2 -0.03-param.priemer/2],'color','black', 'linewidth', 3);
            hold on
        end

        if(param.resultsShow == 1)
            if(modulo == 0)
                del = 1;
                subplot(row,column,order);
                hold on
                for j=1:100:i
                    cm(i) = plot(X(j,param.N+1),X(j,param.N+2),'.','color','red');
                    hold on
                end
                set(gca,'fontsize',18);
                xlabel('X-axis (m)');
                ylabel('Y-axis (m)');
                str = sprintf("Time: %0.2f (s)",(i/100));
                title(str);
                saveas(fig,sprintf('Fig/FIG_%d.fig',saveFig));
                order = order + 1;
            end 
        end

        axis([-dlzkaPotrubia dlzkaPotrubia -0.8 0.8]);
        axis equal
        pause(.001);

        if((i<length(X(:,1))) & (del == 0))
            for lm=1:param.N
                %delete(p(lm));
                %delete(h(lm));
                %delete(m(lm));
    %            delete(cmFct(lm));delete(cmFcn(lm));
                %delete(silX);delete(silY);
                delete(kryt(lm));
            end
        end
        del = 0;
        grid on;      
        axis equal
        axis([-dlzkaPotrubia dlzkaPotrubia -0.8 0.8]);
    end
  
end