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
global M M_p no nw Siw Sor deltafs_shock_p Sw_shock_p  Retar_Norm 

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

Retar_Norm=0.2; % Normalized retadation factor

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
deltafs_p=fw_p./(Sw+Retar_Norm); % slope of fractional flow

% find_shock_p=dfds_p-deltafs_p; % find shock function
Sw_shock_p=fzero('find_shock_PF_Ads',0.9); % find the shock water saturation
Snw_shock_p=(Sw_shock_p-Siw)/(1-Siw-Sor); % find the shock normalized water saturation
fw_shock_p=1/(1+((1-Snw_shock_p)^no/(Snw_shock_p^nw))/M_p); % find the shock water fraction flow
dfds_shock_p=((fw_shock_p^2)/M_p/(1-Siw-Sor))*(((1-Snw_shock_p)^no)/(Snw_shock_p)^nw)*(no/(1-Snw_shock_p)+nw/(Snw_shock_p)); % derivative of water fractional dlow
deltafs_shock_p=fw_shock_p/(Sw_shock_p+Retar_Norm); % slope of fractional flow
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

saveas(figure(1),'Relative Permeability, Water Flooding.tif')

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

saveas(figure(2),'Water Fractional Flow, Water Flooding.tif')

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

saveas(figure(3),'Normalizd Water Fractional Flow, Water Flooding.tif')

figure(4) % Polymer Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)
hold on
S_shock_p=[-Retar_Norm,Sw_shock_p];
F_shock_p=[0,fw_shock_p];
S_shock_p_end=[-Retar_Norm,ER_BT_p_Sw];
F_shock_p_end=[0,1];
plot(S_shock_p,F_shock_p,'r-',S_shock_p_end,F_shock_p_end,'b--','linewidth',1)
hold off

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(4),'Water Fractional Flow, Water and Polymer Flooding,Ads.tif')

figure(5) % Polymer & Water Flooing Fractional Flow Curve and Tangent Line, Water Sturation, Test
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_shock_p=[-Retar_Norm,Sw_shock_p];
F_shock_p=[0,fw_shock_p];
S_shock_p_end=[-Retar_Norm,ER_BT_p_Sw];
F_shock_p_end=[0,1];
plot(S_shock_p,F_shock_p,'r-',S_shock_p_end,F_shock_p_end,'b--','linewidth',1)
plot(S_shock,F_shock,'b-',S_shock_end,F_shock_end,'r--','linewidth',1)
hold off

legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(5),'Water Fractional Flow and shock, Water and Polymer Flooding,Ads.tif')

figure(6) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
hold off

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

text(Sw_shock_p+0.02,fw_shock_p,'Sw_p_f','color','m','fontsize',12)
text(Sw_shock+0.03,fw_shock,'Sw_w_f','color','b','fontsize',12)

saveas(figure(6),'Polymer shock, Polymer Flooding,Ads.tif')

%% Find the intersection, i.e., the oil bank water saturation

Sw_OB=fzero('find_P_OB_Ads',0.5); % find the oil bank water saturation
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
axis([-Retar_Norm 1 0 1])
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

text(Sw_OB+0.02,0.59,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_p+0.021,fw_shock_p-0.02,'Sw_p_f','color','m','fontsize',12)
text(Scw-0.12,0.97,'Sw_i_w','color','k','fontsize',12)

hold off

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(7),'Polymer shock, Water flooded residual oil,Ads.tif')

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
t=0.2;
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

text(0.27,0.40,'Oil Bank','color','b','Rotation',90,'fontsize',14)
text(X_PFb1(1)-0.12,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFt(1)+0.02,0.35,'Sw_O_B','color','g','fontsize',12)
text(X_PFb1(1),0.70,'Sw_p_f','color','m','fontsize',12)
text(X_PFt(1),0.65,'Sw_i_w','color','k','fontsize',12)

X_PFb_Rev=fliplr(X_PFb);
Sw_PFb_Rev=fliplr(Sw_PFb);

X_All=[X_PFb_Rev,X_PFb1,X_PFt,X_PFb2];
Sw_All=[Sw_PFb_Rev,Sw1,Sw_PFt,Sw2];
Pro1=[X_All' Sw_All'];

save PF_Prof_Ads_0.2.mat Pro1

saveas(figure(8),'Saturation profile, Polymer floodings, Water flooded residual oil,Ads.tif')

%% Calculate the Effluent Oil History

figure(9) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
plot(Sw_shock_p,fw_shock_p,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
hold off

legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

text(Sw_OB+0.02,0.60,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_p+0.02,fw_shock_p-0.03,'Sw_p_f','color','m','fontsize',12)
text(Sw_shock+0.01,fw_shock,'Sw_w_f','color','b','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(9),'Polymer shock, compare saturation,Ads.tif')


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

save PFEOR_His_Ads_0.2.mat EOR1

saveas(figure(10),'Polymer EOR,Ads.tif')

%ER_PF4=(Scw-Sw_OB)*(t_OBf3-t_OBf2);
%ER_PF3=[ER_PF1,ER_PF5,ER_PF4];
%figure(11) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
%plot(t_OBf,ER_PF3,'b',t_OBb1,ER_PF2,'b','linewidth',2)
%axis([0 5 0 1-Scw])
%xlabel('Dimensionless Time','fontsize',16)
%ylabel('Recovery Factor','fontsize',16)

%% Compare the effect of adsorption
figure(11) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
load PFEOR_His.mat
t_OBf_All_NoAds=EOR1(:,1);
EOR_All_NoAds=EOR1(:,2);
plot(t_OBf_All_NoAds,EOR_All_NoAds,'b','linewidth',2)

hold on

load PFEOR_His_Ads_0.2.mat
t_OBf_All_Ads=EOR1(:,1);
EOR_All_Ads=EOR1(:,2);
plot(t_OBf_All_Ads,EOR_All_Ads,'g','linewidth',2)

load PFEOR_His_Ads_0.5.mat
t_OBf_All_Ads1=EOR1(:,1);
EOR_All_Ads1=EOR1(:,2);
plot(t_OBf_All_Ads1,EOR_All_Ads1,'r','linewidth',2)

load PFEOR_His_Ads_1.0.mat
t_OBf_All_Ads2=EOR1(:,1);
EOR_All_Ads2=EOR1(:,2);
plot(t_OBf_All_Ads2,EOR_All_Ads2,'k','linewidth',2)

hold off
axis([0 5 0 1-Scw])
xlabel('Dimensionless Time','fontsize',16)
ylabel('Recovery Factor','fontsize',16)

legend({'Polymer Flooding, D_i=0','Polymer Flooding, D_i=0.2','Polymer Flooding, D_i=0.5','Polymer Flooding, D_i=1.0'},'location','NorthWest','fontsize',11)
saveas(figure(11),'Comparison Polymer EOR, all.tif')

figure(12) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
load PF_Prof.mat
X_All=Pro1(:,1);
Sw_All=Pro1(:,2);
plot(X_All,Sw_All,'b','linewidth',2)

hold on

load PF_Prof_Ads_0.2.mat
X_All1=Pro1(:,1);
Sw_All1=Pro1(:,2);
plot(X_All1,Sw_All1,'g','linewidth',2)

load PF_Prof_Ads_0.5.mat
X_All2=Pro1(:,1);
Sw_All2=Pro1(:,2);
plot(X_All2,Sw_All2,'r','linewidth',2)

load PF_Prof_Ads_1.0.mat
X_All3=Pro1(:,1);
Sw_All3=Pro1(:,2);
plot(X_All3,Sw_All3,'k','linewidth',2)

hold off

axis([0 1 0 1])
xlabel('Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)

legend({'Polymer Flooding, D_i=0','Polymer Flooding, D_i=0.2','Polymer Flooding, D_i=0.5','Polymer Flooding, D_i=1.0'},'location','NorthEast','fontsize',11)
saveas(figure(12),'Comparison Saturation Profile of Polymer EOR, all.tif')

