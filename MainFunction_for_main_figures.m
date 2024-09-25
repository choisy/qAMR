function [y_out]=MainFunction()

%{
FIXED PARAMETERS
R00=the R0 of S without drug (R00>1)
pm= growth rate at x= -\infty
p0= growth rate at x=0
p1<p0<pm
%}
pm=10;
p0=.95*pm;
R00=10;
alpha=1;

%{
VARIABLE PARAMETERS
p1_0= ratio p1/p0 range (0,1)
Ratio_cost_ben= log(Delta)/log(1+theta) >0
R_level= level of resistance, a real parameter
k0_factor= the efficiency of the MIC (at x=0) = k0/(1+theta)^R_level
%}

%mutation variance in the phenotypic space
eps=0.01;

%variaance of the initial bacterial pop
sigma0=0.05;
%size of the initial bacterial pop
size_b0=0.05;

%main variable parameter
name='Fig-Main.jpg'; 
kz=0;% just to specify that we are ploting the main figure with all panels
xMinForPlot=0;xMaxForPlot=1;
T=300;% simulation time
p1_0=.5;
k1_0=.01;
k0=10;

%intermediate parameters (as function of main variables parameters)
mu=p0/R00;

p1=p0*p1_0;
k0_1=1/k1_0;

cb1=(1-p1/pm)^-1;
cb0=(1-p0/pm)^-1;

k0_effi0=mu/cb0;
k0_effi1=mu*k0_1/cb1;


Delta=p0*(pm-p1)/(p1*(pm-p0));
Ratio_cost_ben=log(Delta)/log(k0_1);

Curve_S=h_S(Ratio_cost_ben,mu,pm,p0);
Curve_R=h_R(Ratio_cost_ben,mu,pm,p1,k0_1);

k1=k0*k1_0;

xStar=R_evol_level(k0,k1,p1,p0,pm,mu);

[x,time,BB,TotPopB,p,k,R0,BB0,Top]=MainModel(pm,p0,p1,alpha,mu,k0,k1,...
    T,size_b0,sigma0,eps);



figure
set(gcf,'position',[100,100,1300,700])
%axes('fontsize',15)

%FIGURE OF R0  
if kz==0  
    ax=subplot(2,3,1);
else
    ax=subplot(2,2,3);
end
ax.FontSize = 15;
hold on
plot(x,p/mu,'-.m','LineWidth',2)
plot(x,R0,'k','LineWidth',2)
plot([x(1) x(end)], [max(R0) max(R0)], '--k')
plot([xStar xStar], [0 max(R0)],'--k','LineWidth',1)
hold off
xlim([x(1) x(end)])
ax.TickLabelInterpreter = 'latex';
ylabel('Basic reproduction number','Interpreter','latex','fontsize',15)
legend('$R_0^0$','$R_0$','Interpreter','latex','location','northeast')
legend boxoff    
text(xStar,1,'$x^*$','Interpreter','latex','fontsize',15)
if kz==0
    title('A', 'FontSize', 20);
    xlabel('Level of resistance ($x$)','Interpreter','latex','fontsize',15)
else
    xlabel('Level of resistance ($x$)','Interpreter','latex','fontsize',15)
end



%FIGURE initial bacteria and drug efficiency 
if kz==0  
    ax=subplot(2,3,4);
else
    ax=subplot(2,2,4);
end
ax.FontSize = 15;
hold on
plot(x,k,'color',[0.9100    0.4100    0.1700],'LineWidth',2)
Id= x<=1 & x>=-0.15;%Takes from xlim below
area(x,.4*max(k(Id))*BB0/max(BB0),'FaceColor','g','FaceAlpha',.3,'EdgeAlpha',.3)
hold off
text(0.05,.25*max(.4*max(k(Id))*BB0/max(BB0)),'$\longleftarrow$ Initial bacteria','Interpreter','latex','fontsize',12)
xlim([-0.15 1])
xlabel('Level of resistance ($x$)','Interpreter','latex','fontsize',15)
ylabel('Drug efficiency','Interpreter','latex','fontsize',15)
if kz==0
    title('D', 'FontSize', 20);
end


%Plot total bacteria B(t)
if kz==0  
     ax=subplot(2,3,2);
     ax.FontSize = 15;
     plot(log10(time),log10(TotPopB),'LineWidth',2)
     if Top>0
        hold on
        plot([log10(Top) log10(Top)], [log10(min(TotPopB)) log10(max(TotPopB))],'--k','LineWidth',1)
        hold off
     end
     xlabel('$Log_{10}$ of time ($t$)','Interpreter','latex','fontsize',15),
     ylabel('$Log_{10}$ total bacteria','Interpreter','latex','fontsize',15),
     xlim([min(log10(time)) max(log10(time))])
     ylim([log10(min(TotPopB)) log10(max(TotPopB))])
     title('B', 'FontSize', 20);
end 


 %PLot b(t,x)
if kz==0  
    ax=subplot(2,3,[5,6]);
else
    ax=subplot(2,2,[1,2]);
end

ax.FontSize = 15;
colormap parula %autumn %parula
if kz==0  
    VecTime=log10(time);
else
    VecTime=time;
end
surf(VecTime,x,BB,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
set(gca, 'YDir','reverse')
axis tight
view(-136,60)
%colorbar
camlight left
ylabel('$x$','Interpreter','latex','fontsize',15),
zlabel('Bacteria density $b(t,x)$','Interpreter','latex','fontsize',15)
xlim([VecTime(1) VecTime(end)])
if kz==0
    ylim([x(1) xMaxForPlot])
    xlabel('$Log_{10}(t)$','Interpreter','latex','fontsize',15), 
    title('E', 'FontSize', 20);
else
    ylim([xMinForPlot xMaxForPlot])
    xlabel('$t$','Interpreter','latex','fontsize',15), 
    tilName=['Zone ',ZoneID];
    title(tilName, 'FontSize', 20);
end


%FIGURE mean resistance
MeanResis=time;
for n=1:length(time)
    MeanResis(n)=sum(x.*BB(:,n))/TotPopB(n);
end   

if kz==0
    ax=subplot(2,3,3);
    ax.FontSize = 15;
    plot(log10(time),MeanResis,'LineWidth',2)
    xlim([min(log10(time)) max(log10(time))])
    ylim([0,1])
    xlabel('$Log_{10}$ of time ($t$)','Interpreter','latex','fontsize',15)
    ylabel('Average resistance level','Interpreter','latex','fontsize',15)
    ax.TickLabelInterpreter = 'latex';
    title('C', 'FontSize', 20);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function[y]=h_S(z,mu,pm,p0)
       y=mu*(((1-p0/pm)^-1)-z)^-1;
    end

    function[y]=h_R(z,mu,pm,p1,k0_1)
       y=mu*k0_1*(((1-p1/pm)^-1)-z)^-1; 
    end

    function[x,time,BB,TotPopB,p,k,R0,BB0,Top]=MainModel(pm,p0,p1,alpha,mu,k0,k1,...
            T,size_b0,sigma0,eps)
        
        % phenotypic space [0 1]
        xmin=-.5;
        xmax=1.5;  
        % phenotypic space discretization
        Nx=10^2;
        dx=(xmax-xmin)/Nx;
        x=linspace(xmin,xmax,Nx+1)';
        
        %%%% Mutation kernel  
        J=normpdf(dx*(1-Nx:Nx-1),0,eps)';
        
        %functional parameters p and k
        k=k0*(k1/k0).^x;
        pNum=1+(pm/p0-1)*(p0*(pm-p1)/(p1*(pm-p0))).^x;
        p=pm./pNum;
        
        %the R0
        R0=p./(mu+k);
        
        
        % Time discretization
        Nt=10*T;dt=T/Nt;
        time=0:dt:T;
        
        %%%% Initial conditions
        BB0=size_b0* normpdf(x,0,sigma0);
        bb  = BB0;
        
       %%% Time loop (Implicit-explicit Euler scheme)
        BB=zeros(Nx+1,Nt+1);
        TotPopB=zeros(1,Nt+1);
        
        BB(:,1)=bb;
        TotPopB(1)=sum(bb);

%         Solving the system
        for tn=1:Nt 
            GrowthRegulator=(1+TotPopB(tn))^(-alpha);
            IntB=p.*bb;
            LI=GrowthRegulator*dx*(conv(IntB,J,'same'));
            bb=(bb+dt*LI)./(1+dt*mu+dt*k);
            BB(:,tn+1)=bb;
            TotPopB(tn+1)=sum(bb);
        end 
        
        
         %COMPUTING T_op
         Seuil_Top=10^-10;
         if sum(TotPopB)==0
                Top=-3;
         else
                [Row_Top,Col_Top]=find(TotPopB<=Seuil_Top);
                if isempty(Col_Top)
                    Top=-4;
                else
                    id_Top=Col_Top(1);
                    Top=time(id_Top);
                end
         end    
    end

    function[y]=R_evol_level(k0,k1,p1,p0,pm,mu)
       d=k0/k1;b=pm/p0-1;
       a=p0*(pm-p1)/(p1*(pm-p0));
       fun=@(x)(k0*(k1/k0)^x)*log(d)/(mu+(k0*(k1/k0)^x))-...
           b*log(a)*a^x/(1+b*a^x);
       y=fzero(fun,0);
    end
    
end
  