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

global M no nw Siw Sor M_SF no_SF nw_SF Siw_SF Sor_SF deltafs_shock_SF dfds_SF_slug B

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
deltafs_SF=fw_SF./Sw; % slope of fractional flow

% find_shock_p=dfds_p-deltafs_p; % find shock function
Sw_shock_SF=fzero('find_shock_SF',0.9); % find the shock water saturation
Snw_shock_SF=(Sw_shock_SF-Siw_SF)/(1-Siw_SF-Sor_SF); % find the shock normalized water saturation
fw_shock_SF=1/(1+((1-Snw_shock_SF)^no_SF/(Snw_shock_SF^nw_SF))/M_SF); % find the shock water fraction flow
dfds_shock_SF=((fw_shock_SF^2)/M_SF/(1-Siw_SF-Sor_SF))*(((1-Snw_shock_SF)^no_SF)/(Snw_shock_SF)^nw_SF)*(no_SF/(1-Snw_shock_SF)+nw_SF/(Snw_shock_SF)); % derivative of water fractional dlow
deltafs_shock_SF=fw_shock_SF/(Sw_shock_SF); % slope of fractional flow
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

saveas(figure(1),'Relative Permeability, Water Flooding and Surfactant Flooding.tif')

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
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Surfactant Flooding'},'location','northwest','fontsize',12)
hold on
S_shock_SF=[0,Sw_shock_SF];
F_shock_SF=[0,fw_shock_SF];
S_shock_SF_end=[0,ER_BT_SF_Sw];
F_shock_SF_end=[0,1];
plot(S_shock_SF,F_shock_SF,'r-',S_shock_SF_end,F_shock_SF_end,'b--','linewidth',1)
hold off

saveas(figure(4),'Water Fractional Flow, Water and Surfactant Flooding.tif')

figure(5) % Polymer & Water Flooing Fractional Flow Curve and Tangent Line, Water Sturation, Test
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Surfactant Flooding'},'location','northwest','fontsize',12)
hold on
S_shock_SF=[0,Sw_shock_SF];
F_shock_SF=[0,fw_shock_SF];
S_shock_SF_end=[0,ER_BT_SF_Sw];
F_shock_SF_end=[0,1];
plot(S_shock_SF,F_shock_SF,'r-',S_shock_SF_end,F_shock_SF_end,'b--','linewidth',1)
plot(S_shock,F_shock,'b-',S_shock_end,F_shock_end,'r--','linewidth',1)
hold off

saveas(figure(5),'Water Fractional Flow and shock, Water and surfactant Flooding.tif')

figure(6) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line, WF','Tangent Line, SF'},'location','northwest','fontsize',12)

saveas(figure(6),'Surfactant shock, surfactant Flooding.tif')

%% Find the intersection, i.e., the oil bank water saturation

Sw_OB=fzero('find_S_OB',0.5); % find the oil bank water saturation
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
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
S_OB=[Sw_OB,Scw];
F_OB=[fw_OBb,fw_OBf];
plot(S_shock_SF_end,F_shock_SF_end,'r--',S_OB,F_OB,'b--','linewidth',1)


plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Scw,fw_OBf,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

legend({'Water Flooding','Surfactant Flooding','Tangent Line, SF','Shock, Water Front'},'location','northwest','fontsize',12)
% text(2*t,0.6,'Oil Bank','color','b','fontsize',16)

text(Sw_OB-0.12,0.65,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock_SF,fw_shock_SF-0.05,'Sw_s_s','color','m','fontsize',12)
text(Scw-0.12,0.95,'Sw_i_w','color','k','fontsize',12)

hold off

saveas(figure(7),'Surfactant shock, Water flooded residual oil.tif')

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
i_OBt=find(Sw<=(Scw+eps)&Sw>(Scw-eps));
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
t=0.3400;
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

text(2*t-0.05,0.55,'Oil Bank','color','b','fontsize',12)
text(X_SFb1(1)-0.12,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFt(1)+0.02,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFb1(1)+0.02,0.73,'Sw_s_s','color','m','fontsize',12)
text(X_SFt(1),0.65,'Sw_i_w','color','k','fontsize',12)

X_SFb_Rev=fliplr(X_SFb);
Sw_SFb_Rev=fliplr(Sw_SFb);

X_All=[X_SFb_Rev,X_SFb1,X_SFt,X_SFb2];
Sw_All=[Sw_SFb_Rev,Sw1,Sw_SFt,Sw2];
Pro1=[X_All' Sw_All'];

save SF_Prof.mat Pro1

saveas(figure(8),'Saturation profile, Surfactant floodings, Water flooded residual oil,slug.tif')

figure(9) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line compare

plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot(Sw_OB,fw_OBb,'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(Sw_shock,fw_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line, WF','Tangent Line, SF'},'location','northwest','fontsize',12)

text(Sw_OB-0.12,0.65,'Sw_O_B','color','g','fontsize',12)
text(Sw_shock-0.09,fw_shock,'Sw_s','color','c','fontsize',12)
text(Sw_shock_SF,fw_shock_SF-0.05,'Sw_s_s','color','m','fontsize',12)

saveas(figure(9),'Surfactant shock, water shock comparison,slug.tif')

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

save SFEOR_His.mat EOR1

saveas(figure(10),'Surfactant EOR,slug.tif')

%% Slug Injection
t_slug=0.3;
Sw_SF_slug=Sw(i_OBb1+1600);
fw_SF_slug=fw_SF(i_OBb1+1600);
dfds_SF_slug=dfds_SF(i_OBb1+1600);

[A,B]=find_Slug(dfds_SF_slug,Sw_SF_slug,fw_SF_slug);

t_total1=t_slug/A;
X_slug_back=t_slug/B;

figure(11)
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-B 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_SF_slug];
Back_shock_fw=[0,A,fw_SF_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_SF_Slug',0.95); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_SF_slug,fw_SF_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line,S_3','Tangent Line,SF'},'location','NorthWest','fontsize',12)

text(0.03,A,'A','color','k','fontsize',12)
text(-B+0.03,0,'B','color','k','fontsize',12)
text(Sw_shock_SF+0.02,fw_shock_SF-0.05,'S_2','color','m','fontsize',12)
text(Sw_SF_slug+0.02,fw_SF_slug-0.02,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.02,fw_slug_shock,'S_4','color','c','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.17 0.2])

saveas(figure(11),'Surfactant back shock, Surfactant Flooding,Slug,small time.tif')

figure(12)
Front_Oil_Bank_Sw=[Scw,Scw,Sw_OB,Sw_OB];
X_frontOB=t_total1*V_SFf;
X_back_OB=t_total1*dfds_shock_SF;
Front_Oil_Bank_X=[1,t_total1*V_SFf,t_total1*V_SFf,X_back_OB];
plot(Front_Oil_Bank_X,Front_Oil_Bank_Sw,'b','LineWidth',2)
hold on

Surf_Bank_Sw=[Sw_OB,Sw_shock_SF,Sw_SF_slug,Sw_Slug_shock];
Surf_Bank_X=[X_back_OB,X_back_OB,X_slug_back,X_slug_back];
plot(Surf_Bank_X,Surf_Bank_Sw,'b','LineWidth',2)

i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);
dfds_spreading=dfds_SF(i_Spreading_slug:end);
X_Spreading=dfds_spreading.*t_slug;
Sw_spreading=Sw(i_Spreading_slug:end);
X_Spreading_Interp=dfds_spreading.*t_slug*X_slug_back/(X_Spreading(1));
plot(X_Spreading_Interp,Sw_spreading,'b','LineWidth',2)

plot([X_back_OB,X_back_OB],[0,Sw_shock_SF],'r--',[X_slug_back,X_slug_back],[0,Sw_Slug_shock],'r--','LineWidth',2)
% plot([X_Spreading(1),X_slug_back],[Sw_spreading(1),Sw_Slug_shock],'k','LineWidth',2)

plot(X_SFb1(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_SFt(1),Sw(i_OBf),'co','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(X_SFb1(1),Sw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(X_SFt(1),Scw,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(X_slug_back,Sw_SF_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(X_slug_back,Sw_Slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)


text(X_SFb1(1)-0.12,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFt(1)+0.02,0.5,'Sw_O_B','color','g','fontsize',12)
text(X_SFb1(1)+0.02,0.70,'S_2','color','m','fontsize',12)
text(X_SFt(1)-0.05,0.65,'Sw_i_w','color','k','fontsize',12)
text(X_slug_back-0.05,Sw_SF_slug-0.02,'S_3','color','b','fontsize',12)
text(X_slug_back+0.02,Sw_Slug_shock+0.02,'S_4','color','c','fontsize',12)

hold off
axis([0 1 0 1])
xlabel('X_D,Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)
text(0.6,0.55,'Oil Bank','color','b','fontsize',14)
text(0.09,0.30,'Surfactant Bank','color','r','fontsize',14)
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

saveas(figure(12),'Profile Surfactant back shock and oil bank,Slug,small time.tif')

figure(13)
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0.6 0.9 0.8 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_SF_slug];
Back_shock_fw=[0,A,fw_SF_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_SF_Slug',0.95); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_SF_slug,fw_SF_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off

legend({'Water Flooding','Surfactant Flooding','Tangent Line,S_3','Tangent Line,SF'},'location','SouthEast','fontsize',12)

text(Sw_shock_SF+0.005,fw_shock_SF,'S_2','color','m','fontsize',12)
text(Sw_SF_slug+0.005,fw_SF_slug-0.004,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.005,fw_slug_shock,'S_4','color','c','fontsize',12)


saveas(figure(13),'Surfactant back shock, Surfactant Flooding,Slug, magnified,small time.tif')

% Slug Injection
t_slug=0.3;
Sw_slug=Sw(i_OBb1+100);
fw_slug=fw_SF(i_OBb1+100);
dfds_SF_slug=dfds_SF(i_OBb1+100);

[A,B]=find_Slug(dfds_SF_slug,Sw_slug,fw_slug);

t_total=t_slug/A;
X_slug_back=t_slug/B;

figure(14)
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([-B 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_SF_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0.0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0.0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off
legend({'Water Flooding','Surfactant Flooding','Tangent Line,S_3','Tangent Line,SF'},'location','SouthEast','fontsize',12)

text(0.02,A,'A','color','k','fontsize',12)
text(-B+0.02,0.02,'B','color','k','fontsize',12)
text(Sw_shock_SF+0.02,fw_shock_SF-0.05,'S_2','color','m','fontsize',12)
text(Sw_slug+0.02,fw_slug-0.02,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.02,fw_slug_shock,'S_4','color','c','fontsize',12)

ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
hyLabel=get(gca,'YLabel');
% set(hyLabel,'Rotation',90,'VerticalAlignment','bottom')
set(hyLabel,'Rotation',90,'Position',[-0.14 0.2])

saveas(figure(14),'Surfactant back shock, Surfactant Flooding,Slug.tif')


figure(15)
Front_Oil_Bank_Sw=[Scw,Scw,Sw_OB,Sw_OB];
X_frontOB=t_total*V_SFf;
X_back_OB=t_total*dfds_shock_SF;
Front_Oil_Bank_X=[1,t_total*V_SFf,t_total*V_SFf,X_back_OB];
plot(Front_Oil_Bank_X,Front_Oil_Bank_Sw,'b','LineWidth',2)
hold on

Surf_Bank_Sw=[Sw_OB,Sw_shock_SF,Sw_slug,Sw_Slug_shock];
Surf_Bank_X=[X_back_OB,X_back_OB,X_slug_back,X_slug_back];
plot(Surf_Bank_X,Surf_Bank_Sw,'b','LineWidth',2)

i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);
dfds_spreading=dfds_SF(i_Spreading_slug:end);
X_Spreading=dfds_spreading.*t_slug;
Sw_spreading=Sw(i_Spreading_slug:end);
X_Spreading_Interp=dfds_spreading.*t_slug*X_slug_back/(X_Spreading(1));
plot(X_Spreading_Interp,Sw_spreading,'b','LineWidth',2)

plot([X_back_OB,X_back_OB],[0,Sw_shock_SF],'r--',[X_slug_back,X_slug_back],[0,Sw_Slug_shock],'r--','LineWidth',2)
% plot([X_Spreading(1),X_slug_back],[Sw_spreading(1),Sw_Slug_shock],'k','LineWidth',2)
hold off
axis([0 1 0 1])
xlabel('X_D,Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)
text(0.7,0.40,'Oil Bank','color','b','Rotation',90,'fontsize',14)
text(0.25,0.45,'Surfactant Bank','color','r','fontsize',14)
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)

saveas(figure(15),'Profile Surfactant back shock and oil bank,Slug.tif')

figure(16)
plot(Sw,fw,'b',Sw,fw_SF,'r','linewidth',2)
axis([0.6 0.9 0.8 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

Back_shock_Sw=[-B,0,Sw_slug];
Back_shock_fw=[0,A,fw_slug];

hold on
Sw_Slug_shock=fzero('find_S_OB_SF_Slug',0.7); % find the shock normalized water saturation
i_Spreading_slug=find(Sw<=(Sw_Slug_shock+1/2*dSw)&Sw>=(Sw_Slug_shock-1/2*dSw));
fw_slug_shock=fw(i_Spreading_slug);

plot(Back_shock_Sw,Back_shock_fw,'b--',S_shock_SF_end,F_shock_SF_end,'r--','linewidth',1)
plot([Back_shock_Sw(end),Sw_Slug_shock],[Back_shock_fw(end),fw_slug_shock],'b--','linewidth',1)
plot(Sw_shock_SF,fw_shock_SF,'co','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(0,A,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(-B,0,'co','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(Sw_slug,fw_slug,'co','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
plot(Sw_Slug_shock,fw_slug_shock,'co','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',6)
hold off
% legend({'Water Flooding','Polymer Flooding','Tangent Line, WF','Tangent Line, PF'},'location','SouthEast','fontsize',12)
legend({'Water Flooding','Surfactant Flooding','Tangent Line,S_3','Tangent Line,SF'},'location','SouthEast','fontsize',12)
text(Sw_shock_SF+0.005,fw_shock_SF,'S_2','color','m','fontsize',12)
text(Sw_slug+0.005,fw_slug,'S_3','color','b','fontsize',12)
text(Sw_Slug_shock+0.005,fw_slug_shock-0.008,'S_4','color','c','fontsize',12)


saveas(figure(16),'Surfactant back shock, Surfactant Flooding,Slug, magnified.tif')
