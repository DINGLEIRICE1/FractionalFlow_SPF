%% This is the Matlab script for Polymer EOR

%  Analtical Solution; using Buckley-Leverett fractional theory
%  Without considering the interfernce of waves
%  The true residual oil saturation and water/oil relative permeability curves are the same
%  for water flooding and polymer flooding
%  No.5 experiment in Table 5 in SPEJ-179683

%% Define the Variables
clear all
clc
global M M_p no nw Siw Sor deltafs_shock_p Sw_shock_p

Sw=eps:0.01:1; % water saturation
Siw=0.32; % connate water saturation
Sor=0.24; % residual oil saturation
Krw0=0.14; % end point water relative permeability
Kro0=1.0; % end point oil relative permeability
nw=4; % corey component for water
no=2; % corey component for oil
uw=0.5; % water viscosity
up=12; % polymer viscosity
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
fo_p=1-fw_p; % oil fractional flow
dfds_p=((fw_p.^2)./M_p).*(((1-Snw).^no)./(Snw).^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
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


%% Test and Result
%  plot the relative permeability curve

figure(1) % Relative Permeability Curve
plot(Sw(33:77),Krw(33:77),'b',Sw(33:77),Kro(33:77),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Relative Permeability','fontsize',16)
legend({'Krw','Kro'},'fontsize',12)

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

figure(6) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)
hold on
plot(S_shock_end,F_shock_end,'b--',S_shock_p_end,F_shock_p_end,'r--','linewidth',1)
hold off

%% Find the intersection, i.e., the oil bank water saturation

Sw_OB=fzero('find_S_OB',0.5); % find the oil bank water saturation
Snw_OB=(Sw_OB-Siw)/(1-Siw-Sor); % find the oil bank normalized water saturation
fw_OBb=1/(1+((1-Snw_OB)^no/((Snw_OB)^nw))/M); % find the fractional flow of oil bank back

%% Initial Condition
Sio=0.40;
Scw=1-Sio; % connate water saturation
Sncw=(Scw-Siw)/(1-Siw-Sor);
fw_OBf=1/(1+((1-Sncw)^no/((Sncw)^nw))/M); % find the fractional flow of oil bank front

% Test location
figure(7) % Polymer & Water Flooding Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw,'b',Sw,fw_p,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
legend({'Water Flooding','Polymer Flooding'},'location','northwest','fontsize',12)
hold on
S_OB=[Sw_OB,Scw];
F_OB=[fw_OBb,fw_OBf];
plot(S_shock_p_end,F_shock_p_end,'r--',S_OB,F_OB,'b--','linewidth',1)
hold off

%% Interference of Waves
%  none

%% Calculate the Oil Saturation Profile
V_PF=dfds_p;
[deltaf_p_OBb,i_OBb]=max(deltafs_p);
i_OBf=find(Sw>(Sw_OB)&Sw<(Sw_OB+0.01));

% Oil Bank Back
for i=1:length(Sw)
    if i>=i_OBb
        V_PF(i)=dfds_p(i);
    elseif i>=i_OBf
        V_PF(i)=dfds_shock_p;
    else
        V_PF(i)=dfds_p(i);
    end
end
V_PFb=V_PF(i_OBf:end);
Sw_PFb=Sw(i_OBf:end);

% Oil bank Front
i_OBt=find(Sw<=(Scw+eps)&Sw>=(Scw-eps));
V_PFf=(fw_OBf-fw_OBb)/(Scw-Sw(i_OBf));
for i=1:length(Scw)
    if i>=i_OBf && i<=i_OBt
        V_PF(i)=V_PFf;
    else
        V_PF(i)=dfds_p(i);
    end
end
V_PFf=V_PF(i_OBf:i_OBt);
Sw_PFf=Sw(i_OBf:i_OBt);
% Profile
t=0.3;
X_PFb=V_PFb*0.3;
X_PFf=V_PFf*0.3;

figure(8)
plot(X_PFb,Sw_PFb,'b',X_PFf,Sw_PFf,'b','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance','fontsize',16)
ylabel('Water Saturation','fontsize',16)

hold on
X_PFb1=[X_PFb(1),X_PFf(1)];
Sw1=[Sw(i_OBf),Sw(i_OBf)];
X_PFb2=[X_PFf(1),1];
Sw2=[Sw(i_OBt),Sw(i_OBt)];
plot(X_PFb1,Sw1,'b',X_PFb2,Sw2,'b','linewidth',2)
hold off

%% Calculate the Effluent Oil History








