%% This is the Matlab script for Polymer EOR

%  Analtical Solution; using Buckley-Leverett fractional theory
%  Without considering the interfernce of waves
%  The true residual oil saturation and water/oil relative permeability curves are the same
%  for water flooding and polymer flooding
%  No.5 experiment in Table 5 in SPEJ-179683

%% Define the Variables
close all;
clear all;
clc
global M M_p no nw Siw Sor deltafs_shock_p Sw_shock_p dfds_p_slug B

Sw=eps:0.0001:1; % water saturation
dSw=1/(length(Sw)-1); % dSw
Siw=0.15; % connate water saturation
Sor=0.24; % residual oil saturation
Krw0=0.14; % end point water relative permeability
Kro0=0.4; % end point oil relative permeability
nw=4; % corey component for water
no=2; % corey component for oil
uw=0.5; % water viscosity
up=32; % polymer viscosity
uo=72; % oil viscosity
M=Krw0*uo/Kro0/uw; % Mobility ratio

% Normalized water and oil saturation
Snw=Sw; % Normalized water saturation
Sno=Snw; % Normalized oil saturation
Snw=(Sw-Siw)./(1-Siw-Sor); % calculte Normalized water saturation

for i=1:length(Snw)
   if Snw(i)<=0
       Snw(i)=eps;
   elseif Snw(i)>=1
       Snw(i)=1-eps;
   end
end

Sno=1-Snw; % calculate normalized oil saturation

% Relative Permeability Curve
Krw=Snw;
Kro=Sno;
Krw=Krw0*Snw.^nw;
Kro=Kro0*Sno.^no;

% Fractional flow
fw=Snw;
fo=fw;
fw=1./(1+Kro.*uw./(Krw.*uo)); % water fractional flow
fo=1-fw; % oil fractional flow

% dfds(:,i)=((f(:,i).^2)/M(i)).*(((1-S(:)+eps).^n2)./(S(:)+eps).^n1).*(n2./(1-S(:)+eps)+n1./(S(:)+eps));
dfds=((fw.^2)./M).*(((1-Snw).^no)./(Snw).^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
deltafs=fw./(Snw); % slope of fractional flow

% find_shock=dfds-deltafs; % find shock function
Snw_shock=fzero('find_shock',0.5); % find the shock normalized water saturation
fw_shock=1/(1+((1-Snw_shock)^no/(Snw_shock^nw))/M); % find the shock water fraction flow
Sw_shock=Snw_shock*(1-Siw-Sor)+Siw; % find the shock water saturation
dfds_shock=((fw_shock^2)/M)*(((1-Snw_shock)^no)/(Snw_shock)^nw)*(no/(1-Snw_shock)+nw/(Snw_shock)); % derivative of water fractional dlow
deltafs_shock=fw_shock/(Snw_shock); % slope of fractional flow
ER_BT_Snw=Snw_shock-(fw_shock-1)/dfds_shock; % find ER at water B.T.
ER_BT_Sw=ER_BT_Snw*(1-Siw-Sor)+Siw; % find Sw corresponding to ER_BT


%% Polymer Flooding Part
M_p=Krw0*uo/Kro0/up; % Polymer Flooding Mobility ratio
fw_p=Snw;
fo_p=fw;
fw_p=1./(1+Kro.*up./(Krw.*uo));
dfds_p=fw_p;

for i=1:length(Sw)
    fw_p(i)=1/(1+((1-Snw(i))^no/((Snw(i))^nw))/M_p); % fw=1./(1+Kro.*uw./(Krw.*uo));
    dfds_p(i)=((fw_p(i)^2)/M_p/(1-Siw-Sor))*(((1-Snw(i))^no)/(Snw(i))^nw)*(no/(1-Snw(i))+nw/(Snw(i))); % derivative of water fractional dlow
end

% fw_p=1/(1+((1-Snwf)^no/((Snwf)^nw))/M_p); % fw=1./(1+Kro.*uw./(Krw.*uo));

fo_p=1-fw_p; % oil fractional flow
% dfds_p=((fw_p.^2)./M_p).*(((1-Snw).^no)./(Snw).^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
deltafs_p=fw_p./Sw; % slope of fractional flow

% find_shock_p=dfds_p-deltafs_p; % find shock function
Sw_shock_p=fzero('find_shock_p',0.9); % find the shock water saturation
Snw_shock_p=(Sw_shock_p-Siw)/(1-Siw-Sor); % find the shock normalized water saturation
fw_shock_p=1/(1+((1-Snw_shock_p)^no/(Snw_shock_p^nw))/M_p); % find the shock water fraction flow
dfds_shock_p=((fw_shock_p^2)/M_p/(1-Siw-Sor))*(((1-Snw_shock_p)^no)/(Snw_shock_p)^nw)*(no/(1-Snw_shock_p)+nw/(Snw_shock_p)); % derivative of water fractional dlow
deltafs_shock_p=fw_shock_p/(Sw_shock_p); % slope of fractional flow
ER_BT_p_Sw=Sw_shock_p-(fw_shock_p-1)/dfds_shock_p; % find ER at water B.T.
ER_BT_p_Snw=(ER_BT_p_Sw-Siw)/(1-Siw-Sor); % find Sw corresponding to ER_BT

% find the water saturation of polymer bank 

%fw_p=1/(1+((1-Snwf)^no/((Snwf)^nw))/M_p); % fw=1./(1+Kro.*uw./(Krw.*uo));
%dfds_p=((fw_p^2)/M_p/(1-Siw-Sor))*(((1-Snwf)^no)/(Snwf)^nw)*(no/(1-Snwf)+nw/(Snwf)); % derivative of water fractional dlow
%deltafs_p=fw_p/S; % slope of fractional flow
%% Test and Result
%  plot the relative permeability curve

figure(1) % Relative Permeability Curve
plot(Sw(321:761),Krw(321:761),'b',Sw(321:761),Kro(321:761),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Relative Permeability','fontsize',16)
legend({'Krw','Kro'},'fontsize',12)

saveas(figure(1),'Relative Permeability, Water Flooding,Slug.tif')

figure(2) % Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
hold on
S_shock=[Siw,Sw_shock];
F_shock=[0,fw_shock];
S_shock_end=[Siw,ER_BT_Sw];
F_shock_end=[0,1];
plot(S_shock,F_shock,'b-',S_shock_end,F_shock_end,'r--','linewidth',1)
hold off

saveas(figure(2),'Water Fractional Flow, Water Flooding,Slug.tif')

figure(3) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(Snw,fw,'b','linewidth',2)
axis([0 1 0 1])
xlabel('Snw, Normalized Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
hold on
Sn_shock=[0,Snw_shock];
F_shock=[0,fw_shock];
Sn_shock_end=[0,ER_BT_Snw];
F_shock_end=[0,1];
plot(Sn_shock,F_shock,'b-',Sn_shock_end,F_shock_end,'r--','linewidth',1)
hold off

saveas(figure(3),'Normalizd Water Fractional Flow, Water Flooding,Slug.tif')

figure(4) % Polymer Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)
hold on
S_shock_p=[0,Sw_shock_p];
F_shock_p=[0,fw_shock_p];
S_shock_p_end=[0,ER_BT_p_Sw];
F_shock_p_end=[0,1];
plot(S_shock_p,F_shock_p,'r-',S_shock_p_end,F_shock_p_end,'b--','linewidth',1)
hold off

saveas(figure(4),'Water Fractional Flow, Water and Surfactant Flooding,Slug.tif')

figure(5) % Polymer & Water Flooing Fractional Flow Curve and Tangent Line, Water Sturation, Test
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)
hold on
S_shock_p=[0,Sw_shock_p];
F_shock_p=[0,fw_shock_p];
S_shock_p_end=[0,ER_BT_p_Sw];
F_shock_p_end=[0,1];
plot(S_shock_p,F_shock_p,'r-',S_shock_p_end,F_shock_p_end,'b--','linewidth',1)
plot(S_shock,F_shock,'b-',S_shock_end,F_shock_end,'r--','linewidth',1)
hold off

saveas(figure(5),'Water Fractional Flow and shock, Water and Polymer Flooding,Slug.tif')

figure(6) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
hold off
legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)



text(Sw_shock_p+0.02,fw_shock_p,'Sw_p_f','color','m','fontsize',12)
text(Sw_shock+0.03,fw_shock,'Sw_w_f','color','b','fontsize',12)

saveas(figure(6),'Polymer shock, Polymer Flooding,Slug.tif')

%% Find the intersection, i.e., the oil bank water saturation

Sw_OB=fzero('find_P_OB',0.5); % find the oil bank water saturation
Snw_OB=(Sw_OB-Siw)/(1-Siw-Sor); % find the oil bank normalized water saturation
fw_OBb=1/(1+((1-Snw_OB)^no/((Snw_OB)^nw))/M); % find the fractional flow of oil bank back

%% Initial Condition
Sio=0.40; % Initial oil saturation
Scw=1-Sio; % connate water saturation
Sncw=(Scw-Siw)/(1-Siw-Sor);

fw_OBf=1/(1+((1-Sncw)^no/((Sncw)^nw))/M); % find the fractional flow of oil bank front

% Test location
figure(7) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_OB=[Sw_OB,Scw];
F_OB=[fw_OBb,fw_OBf];
plot(S_shock_p_end,F_shock_p_end,'r--',S_OB,F_OB,'b--','linewidth',1)


plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Scw,fw_OBf,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

legend({'Water Flooding','Polymer Flooding','Tangent Line, PF','Shock, Water Front'},'location','SouthEast','fontsize',12)
% text(2*t,0.6,'Oil Bank','color','b','fontsize',16)

text(Sw_OB-0.11,0.50,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_p+0.021,fw_shock_p,'Sw_p_f','color','m','fontsize',12)
text(Scw-0.12,0.97,'Sw_i_w','color','k','fontsize',12)

hold off

saveas(figure(7),'Polymer shock, Water flooded residual oil,Slug.tif')

%% Interference of Waves
%  none

%% Calculate the Oil Saturation Profile
V_PF=dfds_p; %
%[deltaf_p_OBb,i_OBb]=max(deltafs_p); %
i_OBf=find(Sw>=(Sw_OB)&Sw<(Sw_OB+dSw)); %
i_OBb1=find(Sw>=(Sw_shock_p-1/2*dSw)&Sw<=(Sw_shock_p+1/2*dSw));


% Oil Bank Back
for i=1:length(Sw)
    if Sw(i)>=Sw_shock_p
        V_PF(i)=dfds_p(i);
    elseif Sw(i)>=Sw_OB
        V_PF(i)=dfds_shock_p;
    else
        V_PF(i)=dfds_p(i);
    end
end
V_PFb=V_PF(i_OBf:end);
Sw_PFb=Sw(i_OBf:end);

% Oil bank Front
i_OBt=find(Sw<=(Scw+eps)&Sw>=(Scw-eps));
fw_OBt=fw(i_OBt);
fw_OBt1=1/(1+((1-Snw(i_OBt))^no/((Snw(i_OBt))^nw))/M);

V_PFf=(fw_OBf-fw_OBb)/(Scw-Sw(i_OBf));
V_PFt=zeros(1,(i_OBt-i_OBf+1));
V_PFt=V_PFf*ones(1,(i_OBt-i_OBf+1));
Sw_PFt=Sw(i_OBf:i_OBt);

%for i=1:length(Scw)
    %if i>=i_OBf && i<=i_OBt
        %V_PFt(i)=V_PFf;
    %else
        %V_PFt(i)=dfds_p(i);
    %end
%end
%V_PFf=V_PF(i_OBf:i_OBt);
%Sw_PFf=Sw(i_OBf:i_OBt);
% Profile
t=0.4279;
X_PFb=V_PFb*t;
X_PFt=V_PFt*t;

figure(8)
plot(X_PFb,Sw_PFb,'b',X_PFt,Sw_PFt,'b','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)

hold on
X_PFb1=[X_PFb(1),X_PFt(1)];
Sw1=[Sw(i_OBf),Sw(i_OBf)];
X_PFb2=[X_PFt(1),1];
Sw2=[Sw(i_OBt),Sw(i_OBt)];
plot(X_PFb1,Sw1,'b',X_PFb2,Sw2,'b','linewidth',2)

plot(X_PFb1(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_PFt(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_PFb1(1),Sw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(X_PFt(1),Scw,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

hold off

text(0.34,0.40,'Oil Bank','color','b','Rotation',90,'fontsize',14)
text(X_PFb1(1)-0.12,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFt(1)+0.02,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFb1(1),0.70,'Sw_p_f','color','m','fontsize',12)
text(X_PFt(1),0.65,'Sw_i_w','color','k','fontsize',12)

saveas(figure(8),'Saturation profile, Polymer floodings, Water flooded residual oil,Slug.tif')

%% Calculate the Effluent Oil History

figure(9) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
hold off
legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)


text(Sw_OB+0.02,0.45,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_p+0.02,fw_shock_p,'Sw_p_f','color','m','fontsize',12)
text(Sw_shock+0.03,fw_shock,'Sw_w_f','color','b','fontsize',12)

saveas(figure(9),'Polymer shock, compare saturation,Slug.tif')


%% Calculate the Effluent Oil History
t_OBf1=0;
t_OBf2=1/V_PFf; % Breakthrough of oil bank front
t_OBf3=1/dfds_shock_p;


V_PFb1=V_PF(i_OBb1:end);

t_OBb1=V_PFb1;
t_OBb1=1./V_PFb1;

ER_PF1=0;
ER_PF5=0;
ER_PF3=(1-fw_OBb)*(t_OBf3-t_OBf2);
% ER_PF3=(Scw-Sw_OB)*(t_OBf3-t_OBf2);

t_OBf=[t_OBf1,t_OBf2,t_OBf3];
ER_PF=[ER_PF1,ER_PF5,ER_PF3];

ER_PF2=t_OBb1;

for i=1:(length(t_OBb1))
    ER_PF2(i)=Sw(i_OBb1+i-1)-Scw-(fw_p(i_OBb1+i-1)-1)/(dfds_p(i_OBb1+i-1));
end

ER_PF=[ER_PF1,ER_PF5,ER_PF2(1)];

figure(10) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(t_OBf,ER_PF,'b',t_OBb1,ER_PF2,'b','linewidth',2)
axis([0 5 0 1-Scw])
xlabel('Dimensionless Time','fontsize',16)
ylabel('Recovery Factor','fontsize',16)

t_OBf_All=[t_OBf,t_OBb1];
EOR_All=[ER_PF,ER_PF2];
EOR1=[t_OBf_All' EOR_All'];

save PFEOR_His.mat EOR1

saveas(figure(10),'Polymer EOR,Slug.tif')

%% Slug Injection
t_slug=0.3;
Sw_slug=Sw(i_OBb1+500);
fw_slug=fw_p(i_OBb1+500);
dfds_p_slug=dfds_p(i_OBb1+500);

[A,B]=find_Slug(dfds_p_slug,Sw_slug,fw_slug);

t_total=t_slug/A;
X_slug_back=t_slug/B;

figure(11)
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-B 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Polymer Flooding','Tangent Line,S_3','Tangent Line,PF'},'location','NorthWest','fontsize',12)

text(0,A,'A','color','k','fontsize',12)
text(-B,0,'B','color','k','fontsize',12)
text(Sw_shock_p+0.02,fw_shock_p-0.05,'S_2','color','m','fontsize',12)
text(Sw_slug+0.02,fw_slug-0.02,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.02,fw_slug_shock,'S_4','color','c','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])

saveas(figure(11),'Polymer back shock, Polymer Flooding,Slug,small time.tif')

figure(12)
Front_Oil_Bank_Sw=[Scw,Scw,Sw_OB,Sw_OB];
X_frontOB=t_total*V_PFf;
X_back_OB=t_total*dfds_shock_p;
Front_Oil_Bank_X=[1,t_total*V_PFf,t_total*V_PFf,X_back_OB];
plot(Front_Oil_Bank_X,Front_Oil_Bank_Sw,'b','LineWidth',2)
hold on

Poly_Bank_Sw=[Sw_OB,Sw_shock_p,Sw_slug,Sw_Slug_shock];
Poly_Bank_X=[X_back_OB,X_back_OB,X_slug_back,X_slug_back];
plot(Poly_Bank_X,Poly_Bank_Sw,'b','LineWidth',2)

i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);
dfds_spreading=dfds_p(i_Spreading_slug:end);
X_Spreading=dfds_spreading.*t_slug;
Sw_spreading=Sw(i_Spreading_slug:end);
X_Spreading_Interp=dfds_spreading.*t_slug*X_slug_back/(X_Spreading(1));
plot(X_Spreading_Interp,Sw_spreading,'b','LineWidth',2)

plot([X_back_OB,X_back_OB],[0,Sw_shock_p],'r--',[X_slug_back,X_slug_back],[0,Sw_Slug_shock],'r--','LineWidth',2)
% plot([X_Spreading(1),X_slug_back],[Sw_spreading(1),Sw_Slug_shock],'k','LineWidth',2)

plot(X_PFb1(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_PFt(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_PFb1(1),Sw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(X_PFt(1),Scw,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(X_slug_back,Sw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(X_slug_back,Sw_Slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)


text(X_PFb1(1)-0.12,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFt(1)+0.02,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFb1(1)+0.02,0.70,'Sw_s_s','color','m','fontsize',12)
text(X_PFt(1),0.65,'Sw_i_w','color','k','fontsize',12)
text(X_slug_back-0.05,Sw_slug-0.02,'S_3','color','b','fontsize',12)
text(X_slug_back+0.02,Sw_Slug_shock+0.02,'S_4','color','c','fontsize',12)

hold off
axis([0 1 0 1])
xlabel('X_D,Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)
text(0.7,0.40,'Oil Bank','color','b','Rotation',90,'fontsize',14)
text(0.25,0.45,'Polymer Bank','color','r','fontsize',14)
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

saveas(figure(12),'Profile Polymer back shock and oil bank,Slug,small time.tif')


figure(13)
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0.6 0.8 0.8 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Polymer Flooding','Tangent Line,S_3','Tangent Line,PF'},'location','SouthEast','fontsize',12)

text(Sw_shock_p+0.005,fw_shock_p,'S_2','color','m','fontsize',12)
text(Sw_slug+0.005,fw_slug-0.004,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.005,fw_slug_shock,'S_4','color','c','fontsize',12)


saveas(figure(13),'Polymer back shock, Polymer Flooding,Slug, magnified,small time.tif')


% Slug Injection
t_slug=0.3;
Sw_slug=Sw(i_OBb1+100);
fw_slug=fw_p(i_OBb1+100);
dfds_p_slug=dfds_p(i_OBb1+100);

[A,B]=find_Slug(dfds_p_slug,Sw_slug,fw_slug);

t_total1=t_slug/A;
X_slug_back=t_slug/B;

figure(14)
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-B 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0.0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0.0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off
legend({'Water Flooding','Polymer Flooding','Tangent Line,S_3','Tangent Line,PF'},'location','SouthEast','fontsize',12)

text(0.02,A,'A','color','k','fontsize',12)
text(-B+0.02,0.02,'B','color','k','fontsize',12)
text(Sw_shock_p+0.02,fw_shock_p-0.05,'S_2','color','m','fontsize',12)
text(Sw_slug+0.02,fw_slug-0.02,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.02,fw_slug_shock,'S_4','color','c','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])

saveas(figure(14),'Polymer back shock, Polymer Flooding,Slug.tif')

figure(15)
Front_Oil_Bank_Sw=[Scw,Scw,Sw_OB,Sw_OB];
X_frontOB=t_total1*V_PFf;
X_back_OB=t_total1*dfds_shock_p;
Front_Oil_Bank_X=[1,t_total1*V_PFf,t_total1*V_PFf,X_back_OB];
plot(Front_Oil_Bank_X,Front_Oil_Bank_Sw,'b','LineWidth',2)
hold on

Poly_Bank_Sw=[Sw_OB,Sw_shock_p,Sw_slug,Sw_Slug_shock];
Poly_Bank_X=[X_back_OB,X_back_OB,X_slug_back,X_slug_back];
plot(Poly_Bank_X,Poly_Bank_Sw,'b','LineWidth',2)

i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);
dfds_spreading=dfds_p(i_Spreading_slug:end);
X_Spreading=dfds_spreading.*t_slug;
Sw_spreading=Sw(i_Spreading_slug:end);
X_Spreading_Interp=dfds_spreading.*t_slug*X_slug_back/(X_Spreading(1));
plot(X_Spreading_Interp,Sw_spreading,'b','LineWidth',2)

plot([X_back_OB,X_back_OB],[0,Sw_shock_p],'r--',[X_slug_back,X_slug_back],[0,Sw_Slug_shock],'r--','LineWidth',2)
% plot([X_Spreading(1),X_slug_back],[Sw_spreading(1),Sw_Slug_shock],'k','LineWidth',2)
hold off
axis([0 1 0 1])
xlabel('X_D,Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)
text(0.7,0.40,'Oil Bank','color','b','Rotation',90,'fontsize',14)
text(0.25,0.45,'Polymer Bank','color','r','fontsize',14)
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

saveas(figure(15),'Profile Polymer back shock and oil bank,Slug.tif')


figure(16)
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0.6 0.8 0.8 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)
legend({'Water Flooding','Polymer Flooding','Tangent Line,S_3','Tangent Line,PF'},'location','SouthEast','fontsize',12)
text(Sw_shock_p+0.005,fw_shock_p,'S_2','color','m','fontsize',12)
text(Sw_slug+0.005,fw_slug,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.005,fw_slug_shock,'S_4','color','c','fontsize',12)


saveas(figure(16),'Polymer back shock, Polymer Flooding,Slug, magnified.tif')
