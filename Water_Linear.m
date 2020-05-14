%% This is the Matlab script for Linear Water Flooding
%  Analtical Solution; using Buckley-Leverett fractional theory
%  Without considering the interfernce of waves
%  The true residual oil saturation and water/oil relative permeability curves are the same
%  Without Considering the gravity, capillary pressure or dispersion

%% Define the Variables
clear all
clc
global M no nw Siw Sor

Sw=0:0.001:1; % water saturation
Siw=0.32; % connate water saturation
Sor=0.24; % residual oil saturation
Krw0=0.14; % end point water relative permeability
Kro0=1.0; % end point oil relative permeability
nw=4; % corey component for water
no=2; % corey component for oil
uw=0.5; % water viscosity
up=12; % polymer viscosity
uo=10; % oil viscosity
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
% slope of fractional flow
dfds=((fw.^2)./M).*(((1-Snw).^no)./(Snw).^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
deltafs=fw./(Snw);
% find_shock=dfds-deltafs; % find shock function
Snw_shock=fzero('find_shock',0.5); % find the shock normalized water saturation
fw_shock=1/(1+((1-Snw_shock)^no/(Snw_shock^nw))/M); % find the shock water fraction flow
Sw_shock=Snw_shock*(1-Siw-Sor)+Siw; % find the shock water saturation
dfds_shock=((fw_shock^2)/M)*(((1-Snw_shock)^no)/(Snw_shock)^nw)*(no/(1-Snw_shock)+nw/(Snw_shock)); % derivative of water fractional dlow
deltafs_shock=fw_shock/(Snw_shock); % slope of fractional flow

% Calculate the Effluent Oil History
t_BT=1/deltafs_shock; % Find B.T. time of water
ER_BT_Snw=Snw_shock-(fw_shock-1)/dfds_shock; % find ER at water B.T.
ER_BT_Sw=ER_BT_Snw*(1-Siw-Sor)+Siw; % find Sw corresponding to ER_BT

ER_WF_Snw=Snw;
Wcut_WF_Snw=Snw;

% Saturation velocity
V_Snw=Snw;

% calcuate the effluent history
i_Sw_shock=find(Sw>=(Sw_shock)&Sw<(Sw_shock+0.001));
for i=1:length(Snw)
    if i<i_Sw_shock
       V_Snw(i)=deltafs_shock;
       Wcut_WF_Snw(i)=0;
    else
       V_Snw(i)=dfds(i);
       Wcut_WF_Snw(i)=fw(i);
    end
end

T_WF_Snw=1./V_Snw;
n_ER_WF_Snw=T_WF_Snw(1)/(i_Sw_shock-1);
ER_WF_Snw1=0:n_ER_WF_Snw:T_WF_Snw(1);

for i=1:length(Snw)
    if i<i_Sw_shock
       ER_WF_Snw(i)=ER_WF_Snw1(i);
       T_WF_Snw(i)=ER_WF_Snw(i);
    else
       ER_WF_Snw(i)=Snw(i)-(fw(i)-1)/dfds(i);
    end
end


%% Test and Result
%  plot the relative permeability curve

figure(1) % Relative Permeability Curve
plot(Sw(320:760),Krw(320:760),'b',Sw(320:760),Kro(320:760),'r','linewidth',2)
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

figure(4) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw,ER_WF_Snw,'b','linewidth',2)
axis([0 5 0 1])
xlabel('Dimentionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)

figure(5) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw,Wcut_WF_Snw,'b','linewidth',2)
axis([0 5 0 1])
xlabel('Dimentionless Time, MPV','fontsize',16)
ylabel('Water Cut','fontsize',16)

