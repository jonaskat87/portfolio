close all
ms = 1 ;
mt = 0.9 ;
mj = 1 ;
G = 0.85 ;
CI = [1 -0.3 1 1 2 0 0 0 0 0 0 0];
to = 0;
tf = 200;
tol = 1e-8;


F = @(t,Y) [Y(7);Y(8);Y(9);Y(10);Y(11);Y(12);G*(mt*(Y(3)-Y(1))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)+mj*(Y(5)-Y(1))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2));G*(mt*(Y(4)-Y(2))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)+mj*(Y(6)-Y(2))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2));G*(ms*(Y(1)-Y(3))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)+mj*(Y(5)-Y(3))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));G*(ms*(Y(2)-Y(4))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)+mj*(Y(6)-Y(4))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));G*(ms*(Y(1)-Y(5))/((Y(1)-Y(5))^2+(Y(6)-Y(2))^2)^(3/2)+mt*(Y(3)-Y(5))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));G*(ms*(Y(2)-Y(6))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2)+mt*(Y(4)-Y(6))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2))];
options = odeset('RelTol', tol, 'AbsTol', tol);
[t,Y]=ode45(F,[to,tf],CI,options);

%plot(Y(:,1),Y(:,2),'r',Y(:,3),Y(:,4),'g',Y(:,5),Y(:,6),'b')
yo = Y(1) ;
dto = 0.3 ;
plotmax = 100 ;
T=to ;

xmin = min([Y(:,1);Y(:,3);Y(:,5)]);
xmax = max([Y(:,1);Y(:,3);Y(:,5)]);
ymin = min([Y(:,2);Y(:,4);Y(:,6)]);
ymax = max([Y(:,2);Y(:,4);Y(:,6)]);

figure('Position',[0 0 1550 800])
for i = 1:length(Y(:,1))
    
    dt = abs(Y(i,1)-yo)/abs(Y(i,7));
    if dt >= dto
        %plot(Y(i,1),Y(i,2),'--or',Y(i,3),Y(i,4),'--og',Y(i,5),Y(i,6),'--ob')
        if i > plotmax
            plot(Y(i-plotmax:i,1),Y(i-plotmax:i,2),'r',Y(i,1),Y(i,2),'-or',Y(i-plotmax:i,3),Y(i-plotmax:i,4),'g',Y(i,3),Y(i,4),'-og',Y(i-plotmax:i,5),Y(i-plotmax:i,6),'b',Y(i,5),Y(i,6),'-ob')
            title(sprintf('t=%f',t(i)))
            T=[T;t(i)];
        else
            
            plot(Y(1:i,1),Y(1:i,2),'r',Y(i,1),Y(i,2),'-or',Y(1:i,3),Y(1:i,4),'g',Y(i,3),Y(i,4),'-og',Y(1:i,5),Y(1:i,6),'b',Y(i,5),Y(i,6),'-ob')
            title(sprintf('t=%f',t(i)))
            T=[T;t(i)];
        end
            
        axis([xmin xmax ymin ymax])
        yo = Y(i,1) ; 
        vo = Y(i,7) ;
    end
    pause(0.001)
end
X=[0:1:length(T)-1];
%figure(2)
%plot(X,T)
%plot(Y(:,1),Y(:,2),'r',Y(:,3),Y(:,4),'g',Y(:,5),Y(:,6),'b')