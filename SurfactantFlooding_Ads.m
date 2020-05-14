%% This is the Matlab script for Surfactant EOR

%  Analtical Solution; using Buckley-Leverett fractional theory
%  Without considering the interfernce of waves
%  The true residual oil saturation and water/oil relative permeability curves are the same
%  for water flooding and Surfactant flooding
%  No.5 experiment in Table 5 in SPEJ-179683

%% Define the Variables for water flooding
close all;
clear all;
clc

global M no nw Siw Sor M_SF no_SF nw_SF Siw_SF Sor_SF deltafs_shock_SF Retar_Norm

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
uo=5; % oil viscosity
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

%% Define the Variables for Surfactant Flooding

Siw_SF=0.15; % connate water saturation
Sor_SF=0.10; % residual oil saturation
Krw0_SF=0.4; % end point water relative permeability
Kro0_SF=1.0; % end point oil relative permeability
nw_SF=4; % corey component for water
no_SF=1.5; % corey component for oil
M_SF=Krw0_SF*uo/Kro0_SF/uw; % Mobility ratio

Retar_Norm=0.2; % Normalized retadation factor

% Normalized water and oil saturation
Snw_SF=Sw; % Normalized water saturation
Sno_SF=Snw_SF; % Normalized oil saturation
Snw_SF=(Sw-Siw_SF)./(1-Siw_SF-Sor_SF); % calculte Normalized water saturation

for i=1:length(Snw_SF)
   if Snw_SF(i)<=0
       Snw_SF(i)=eps;
   elseif Snw_SF(i)>=1
       Snw_SF(i)=1-eps;
   end
end

Sno_SF=1-Snw_SF; % calculate normalized oil saturation

% Relative Permeability Curve
Krw_SF=Snw_SF;
Kro_SF=Sno_SF;
Krw_SF=Krw0_SF*Snw_SF.^nw_SF;
Kro_SF=Kro0_SF*Sno_SF.^no_SF;

% Water Flooding Fractional flow
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


%% Surfactant Flooding Fractional Flow
M_SF=Krw0_SF*uo/Kro0_SF/uw; % Surfactant Flooding Mobility ratio
fw_SF=Snw_SF;
fo_SF=fw_SF;

dfds_SF=fw_SF;

for i=1:length(Sw)
    fw_SF(i)=1/(1+((1-Snw_SF(i))^no_SF/((Snw_SF(i))^nw_SF))/M_SF); % fw=1./(1+Kro.*uw./(Krw.*uo));
    dfds_SF(i)=((fw_SF(i)^2)/M_SF/(1-Siw_SF-Sor_SF))*(((1-Snw_SF(i))^no_SF)/(Snw_SF(i))^nw_SF)*(no_SF/(1-Snw_SF(i))+nw_SF/(Snw_SF(i))); % derivative of water fractional dlow
end

fo_SF=1-fw_SF; % oil fractional flow
% dfds_SF=((fw_SF.^2)./M_SF).*(((1-Snw_SF).^no_SF)./(Snw_SF).^nw_SF).*(no_SF./(1-Snw_SF)+nw_SF./(Snw_SF)); % derivative of water fractional dlow
deltafs_SF=fw_SF./(Sw+Retar_Norm); % slope of fractional flow

% find_shock_p=dfds_p-deltafs_p; % find shock function
Sw_shock_SF=fzero('find_shock_SF_Ads',0.9); % find the shock water saturation
Snw_shock_SF=(Sw_shock_SF-Siw_SF)/(1-Siw_SF-Sor_SF); % find the shock normalized water saturation
fw_shock_SF=1/(1+((1-Snw_shock_SF)^no_SF/(Snw_shock_SF^nw_SF))/M_SF); % find the shock water fraction flow
dfds_shock_SF=((fw_shock_SF^2)/M_SF/(1-Siw_SF-Sor_SF))*(((1-Snw_shock_SF)^no_SF)/(Snw_shock_SF)^nw_SF)*(no_SF/(1-Snw_shock_SF)+nw_SF/(Snw_shock_SF)); % derivative of water fractional dlow
deltafs_shock_SF=fw_shock_SF/(Sw_shock_SF+Retar_Norm); % slope of fractional flow
ER_BT_SF_Sw=Sw_shock_SF-(fw_shock_SF-1)/dfds_shock_SF; % find ER at water B.T.
ER_BT_SF_Snw=(ER_BT_SF_Sw-Siw_SF)/(1-Siw_SF-Sor_SF); % find Sw corresponding to ER_BT

% find the water saturation of polymer bank 


%% Test and Result
%  plot the relative permeability curve

figure(1) % Relative Permeability Curve
plot(Sw,Krw,'b',Sw,Kro,'r','linewidth',2)
hold on
plot(Sw,Krw_SF,'b--',Sw,Kro_SF,'r--','linewidth',2)
hold off
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Relative Permeability','fontsize',16)
legend({'Krw','Kro','Krw,SF','Kro, SF'},'fontsize',12)

saveas(figure(1),'Relative Permeability, Water Flooding and Surfactant Flooding, with Ads.tif')

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

saveas(figure(2),'Water Fractional Flow, Water Flooding, with Ads.tif')

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

saveas(figure(3),'Normalizd Water Fractional Flow, Water Flooding, with Ads.tif')

figure(4) % Polymer Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_shock_SF=[-Retar_Norm,Sw_shock_SF];
F_shock_SF=[0,fw_shock_SF];
S_shock_SF_end=[-Retar_Norm,ER_BT_SF_Sw];
F_shock_SF_end=[0,1];
% plot(S_shock_SF,F_shock_SF,'r-',S_shock_SF_end,F_shock_SF_end,'b--','linewidth',1)
plot(S_shock_SF,F_shock_SF,'r-','linewidth',1)

legend({'Water Flooding','Surfactant Flooding','Tangent Line'},'location','southeast','fontsize',11)

hold off

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(4),'Water Fractional Flow, Water and Surfactant Flooding, with Ads.tif')

figure(5) % Polymer & Water Flooing Fractional Flow Curve and Tangent Line, Water Sturation, Test
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_shock_SF=[-Retar_Norm,Sw_shock_SF];
F_shock_SF=[0,fw_shock_SF];
S_shock_SF_end=[-Retar_Norm,ER_BT_SF_Sw];
F_shock_SF_end=[0,1];
plot(S_shock_SF,F_shock_SF,'r-',S_shock_SF_end,F_shock_SF_end,'b--','linewidth',1)
plot(S_shock,F_shock,'b-',S_shock_end,F_shock_end,'r--','linewidth',1)
hold off

legend({'Water Flooding','Surfactant Flooding'},'location','SouthEast','fontsize',11)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(5),'Water Fractional Flow and shock, Water and surfactant Flooding, with Ads.tif')

figure(6) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line, WF','Tangent Line, SF'},'location','SouthEast','fontsize',11)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(6),'Surfactant shock, surfactant Flooding, with Ads.tif')

%% Find the intersection, i.e., the oil bank water saturation

Sw_OB=fzero('find_S_OB_Ads',0.5); % find the oil bank water saturation
Snw_OB=(Sw_OB-Siw)/(1-Siw-Sor); % find the oil bank normalized water saturation
fw_OBb=1/(1+((1-Snw_OB)^no/((Snw_OB)^nw))/M); % find the fractional flow of oil bank back

%% Initial Condition
Sio=0.40; % Initial oil saturation
Scw=1-Sio; % connate water saturation
Sncw=(Scw-Siw)/(1-Siw-Sor);

fw_OBf=1/(1+((1-Sncw)^no/((Sncw)^nw))/M); % find the fractional flow of oil bank front

% Test location
figure(7) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_OB=[Sw_OB,Scw];
F_OB=[fw_OBb,fw_OBf];
plot(S_shock_SF_end,F_shock_SF_end,'r--',S_OB,F_OB,'b--','linewidth',1)


plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Scw,fw_OBf,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

legend({'Water Flooding','Surfactant Flooding','Tangent Line, SF','Shock, Water Front'},'location','SouthEast','fontsize',11)
% text(2*t,0.6,'Oil Bank','color','b','fontsize',16)

text(Sw_OB-0.135,0.75,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_SF,fw_shock_SF-0.05,'Sw_s_s','color','m','fontsize',12)
text(Scw-0.12,0.95,'Sw_i_w','color','k','fontsize',12)

hold off

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(7),'Surfactant shock, Water flooded residual oil, with Ads.tif')

%% Interference of Waves
%  none

%% Calculate the Oil Saturation Profile
V_SF=dfds_SF; %
%[deltaf_p_OBb,i_OBb]=max(deltafs_p); %
i_OBf=find(Sw>=(Sw_OB)&Sw<(Sw_OB+dSw)); %
i_OBb1=find(Sw>=(Sw_shock_SF-1/2*dSw)&Sw<=(Sw_shock_SF+1/2*dSw));

% Oil Bank Back
for i=1:length(Sw)
    if Sw(i)>=Sw_shock_SF
        V_SF(i)=dfds_SF(i);
    elseif Sw(i)>=Sw_OB
        V_SF(i)=dfds_shock_SF;
    else
        V_SF(i)=dfds_SF(i);
    end
end

V_SFb=V_SF(i_OBf:end);
Sw_SFb=Sw(i_OBf:end);

% Oil bank Front
i_OBt=find(Sw<=(Scw+eps)&Sw>=(Scw-eps));
fw_OBt=fw(i_OBt);
fw_OBt1=1/(1+((1-Snw(i_OBt))^no/((Snw(i_OBt))^nw))/M);

V_SFf=(fw_OBf-fw_OBb)/(Scw-Sw(i_OBf));
V_SFt=zeros(1,(i_OBt-i_OBf+1));
V_SFt=V_SFf*ones(1,(i_OBt-i_OBf+1));
Sw_SFt=Sw(i_OBf:i_OBt);

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
X_SFb=V_SFb*t;
X_SFt=V_SFt*t;

figure(8)
plot(X_SFb,Sw_SFb,'b',X_SFt,Sw_SFt,'b','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)

hold on
X_SFb1=[X_SFb(1),X_SFt(1)];
Sw1=[Sw(i_OBf),Sw(i_OBf)];
X_SFb2=[X_SFt(1),1];
Sw2=[Sw(i_OBt),Sw(i_OBt)];
plot(X_SFb1,Sw1,'b',X_SFb2,Sw2,'b','linewidth',2)

plot(X_SFb1(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_SFt(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_SFb1(1),Sw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(X_SFt(1),Scw,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

hold off

text(2*t-0.10,0.55,'Oil Bank','color','b','fontsize',12)
text(X_SFb1(1)-0.12,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFt(1)+0.02,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFb1(1)+0.02,0.73,'Sw_s_s','color','m','fontsize',12)
text(X_SFt(1),0.65,'Sw_i_w','color','k','fontsize',12)

X_SFb_Rev=fliplr(X_SFb);
Sw_SFb_Rev=fliplr(Sw_SFb);

X_All=[X_SFb_Rev,X_SFb1,X_SFt,X_SFb2];
Sw_All=[Sw_SFb_Rev,Sw1,Sw_SFt,Sw2];
Pro1=[X_All' Sw_All'];

save SF_Prof_Ads_0.2.mat Pro1

saveas(figure(8),'Saturation profile, Surfactant floodings, Water flooded residual oil, with Ads.tif')

%figure (12)
%plot(X_All,Sw_All,'k','linewidth',2)
%axis([0 1 0 1])

figure(9) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line compare

plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-Retar_Norm 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line, WF','Tangent Line, SF'},'location','SouthEast','fontsize',11)

text(Sw_OB-0.13,0.75,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock-0.10,fw_shock,'Sw_s','color','c','fontsize',12)
text(Sw_shock_SF,fw_shock_SF-0.05,'Sw_s_s','color','m','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])
% set(get(gca,'Ylabel'),'Rotation',90)

saveas(figure(9),'Surfactant shock, water shock comparison, with Ads.tif')

%% Calculate the Effluent Oil History
t_OBf1=0;
t_OBf2=1/V_SFf; % Breakthrough of oil bank front
t_OBf3=1/dfds_shock_SF;


V_SFb1=V_SF(i_OBb1:end);

t_OBb1=V_SFb1;
t_OBb1=1./V_SFb1;

ER_SF1=0;
ER_SF5=0;
ER_SF3=(1-fw_OBb)*(t_OBf3-t_OBf2);
% ER_PF3=(Scw-Sw_OB)*(t_OBf3-t_OBf2);

t_OBf=[t_OBf1,t_OBf2,t_OBf3];
ER_SF=[ER_SF1,ER_SF5,ER_SF3];

ER_SF2=t_OBb1;

for i=1:(length(t_OBb1))
    ER_SF2(i)=Sw(i_OBb1+i-1)-Scw-(fw_SF(i_OBb1+i-1)-1)/(dfds_SF(i_OBb1+i-1));
end

ER_SF=[ER_SF1,ER_SF5,ER_SF2(1)];

figure(10) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(t_OBf,ER_SF,'b',t_OBb1,ER_SF2,'b','linewidth',2)
axis([0 5 0 1-Scw])
xlabel('Dimensionless Time','fontsize',16)
ylabel('Recovery Factor','fontsize',16)

t_OBf_All=[t_OBf,t_OBb1];
EOR_All=[ER_SF,ER_SF2];
EOR1=[t_OBf_All' EOR_All'];

save SFEOR_His_Ads_0.2.mat EOR1

saveas(figure(10),'Surfactant EOR, With Ads.tif')

%% Compare the effect of adsorption
figure(11) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
load SFEOR_His.mat
t_OBf_All_NoAds=EOR1(:,1);
EOR_All_NoAds=EOR1(:,2);
plot(t_OBf_All_NoAds,EOR_All_NoAds,'b','linewidth',2)

hold on

load SFEOR_His_Ads_0.2.mat
t_OBf_All_Ads=EOR1(:,1);
EOR_All_Ads=EOR1(:,2);
plot(t_OBf_All_Ads,EOR_All_Ads,'g','linewidth',2)

load SFEOR_His_Ads_0.5.mat
t_OBf_All_Ads1=EOR1(:,1);
EOR_All_Ads1=EOR1(:,2);
plot(t_OBf_All_Ads1,EOR_All_Ads1,'r','linewidth',2)

load SFEOR_His_Ads_1.0.mat
t_OBf_All_Ads2=EOR1(:,1);
EOR_All_Ads2=EOR1(:,2);
plot(t_OBf_All_Ads2,EOR_All_Ads2,'k','linewidth',2)

hold off
axis([0 5 0 1-Scw])
xlabel('Dimensionless Time','fontsize',16)
ylabel('Recovery Factor','fontsize',16)

legend({'Surfactant Flooding, D_i=0','Surfactant Flooding, D_i=0.2','Surfactant Flooding, D_i=0.5','Surfactant Flooding, D_i=1.0'},'location','NorthWest','fontsize',11)
saveas(figure(11),'Comparison Surfactant EOR, all.tif')

figure(12) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
load SF_Prof.mat
X_All=Pro1(:,1);
Sw_All=Pro1(:,2);
plot(X_All,Sw_All,'b','linewidth',2)

hold on

load SF_Prof_Ads_0.2.mat
X_All1=Pro1(:,1);
Sw_All1=Pro1(:,2);
plot(X_All1,Sw_All1,'g','linewidth',2)

load SF_Prof_Ads_0.5.mat
X_All2=Pro1(:,1);
Sw_All2=Pro1(:,2);
plot(X_All2,Sw_All2,'r','linewidth',2)

load SF_Prof_Ads_1.0.mat
X_All3=Pro1(:,1);
Sw_All3=Pro1(:,2);
plot(X_All3,Sw_All3,'k','linewidth',2)

hold off

axis([0 1 0 1])
xlabel('Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)

legend({'Surfactant Flooding, D_i=0','Surfactant Flooding, D_i=0.2','Surfactant Flooding, D_i=0.5','Surfactant Flooding, D_i=1.0'},'location','NorthEast','fontsize',11)
saveas(figure(12),'Comparison Saturation Profile of Surfactant EOR, all.tif')
