%% This is the Matlab script for Linear Water Flooding

%  Bug: for piston displacement, the water cut and oil cut history are not correct 
%  Compare the effect of mobility ratio on oil recovery efficiency
%  Analtical Solution; using Buckley-Leverett fractional theory
%  Without considering the interfernce of waves
%  The true residual oil saturation and water/oil relative permeability curves are the same
%  Without Considering the gravity, capillary pressure or dispersion

%% Define the Variables
clear all
clc
global no nw Siw Sor

Sw=0:0.001:1; % water saturation
dSw=1/(length(Sw)-1);
Siw=0.32; % connate water saturation
Sor=0.24; % residual oil saturation
Krw0=0.14; % end point water relative permeability
Kro0=1.0; % end point oil relative permeability
nw=4; % corey component for water
no=2; % corey component for oil
uw=0.5; % water viscosity
up=12; % polymer viscosity
uo=[0.5 72 1200]; % oil viscosity

for i=1:length(uo)
    M(i)=Krw0*uo(i)/Kro0/uw; % Mobility ratio
end


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

Krw=Krw0*Snw.^nw;
Kro=Kro0*Sno.^no;

% Fractional flow
fw=ones(length(M),length(Krw));
fo=fw;
for i=1:length(uo)
    fw(i,:)=1./(1+Kro.*uw./(Krw.*uo(i))); % water fractional flow
    fo=1-fw; % oil fractional flow
end

% dfds(:,i)=((f(:,i).^2)/M(i)).*(((1-S(:)+eps).^n2)./(S(:)+eps).^n1).*(n2./(1-S(:)+eps)+n1./(S(:)+eps));
% slope of fractional flow
dfds=fw;
deltafs=fw;
for i=1:length(uo)
    dfds(i,:)=((fw(i,:).^2)./M(i)).*(((1-Snw)).^no)./(Snw.^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
    deltafs(i,:)=fw(i,:)./(Snw);
end

% find_shock=dfds-deltafs; % find shock function
Snw_shock=ones(1,length(M));
fw_shock=Snw_shock;
Sw_shock=fw_shock;
dfds_shock=fw_shock;
deltafs_shock=fw_shock;
t_BT=fw_shock;
ER_BT_Snw=fw_shock;

for i=1:length(uo)
    % find_shock=dfds-deltafs; % find shock function
    
    Snw_shock(1,i)=fzero('find_shock_Com',0.7,[],M(i)); % find the shock normalized water saturation
    fw_shock(1,i)=1/(1+((1-Snw_shock(1,i))^no/(Snw_shock(1,i)^nw))/M(i)); % find the shock water fraction flow
    Sw_shock(1,i)=Snw_shock(1,i)*(1-Siw-Sor)+Siw; % find the shock water saturation
    dfds_shock(1,i)=((fw_shock(1,i)^2)/M(i))*(((1-Snw_shock(1,i))^no)/(Snw_shock(1,i))^nw)*(no/(1-Snw_shock(1,i))+nw/(Snw_shock(1,i))); % derivative of water fractional dlow
    deltafs_shock(1,i)=fw_shock(1,i)/(Snw_shock(1,i)); % slope of fractional flow
    % Calculate the Effluent Oil History

    t_BT(1,i)=1/deltafs_shock(1,i);
    ER_BT_Snw(1,i)=Snw_shock(1,i)-(fw_shock(1,i)-1)/dfds_shock(1,i); % find ER at water B.T.
    ER_BT_Sw(1,i)=ER_BT_Snw(1,i)*(1-Siw-Sor)+Siw; % find Sw corresponding to ER_BT
    
    i_Sw_shock(i)=find(Sw>(Sw_shock(1,i))&Sw<(Sw_shock(1,i)+dSw));
    
    if Sw_shock(1,i)>=1-Sor-dSw
        Snw_shock(1,i)=1;
        Sw_shock(1,i)=1-Sor;
        fw_shock(1,i)=1;
        dfds_shock(1,i)=1;
        deltafs_shock(1,i)=1;
        i_Sw_shock(i)=find(Sw>(1-Sor-dSw) & Sw<(1-Sor+dSw));
    end
end

% Calculate the Recovery Factor

ER_WF_Snw=fw;
Wcut_WF_Snw=fw;
Ocut_WF_Snw=fw;
T_WF_Snw=fw;
% Saturation velocity
V_Snw=fw;
% calcuate the effluent history

for i=1:length(uo) 
    for j=1:length(Snw)
        if j<i_Sw_shock(i)
            V_Snw(i,j)=deltafs_shock(1,i);
            Wcut_WF_Snw(i,j)=0;
        else
            V_Snw(i,j)=dfds(i,j);
            Wcut_WF_Snw(i,j)=fw(i,j);
        end
    end
    
    if Snw(1,i)==1
       for j=1:length(Snw)
           if j<i_Sw_shock
               V_Snw(i,j)=1;
               Wcut_WF_Snw(i,j)=0;
           else
               % V_Snw(i,j)=1;
               Wcut_WF_Snw(i,j)=1;
           end
       end
    end

end

Ocut_WF_Snw=1-Wcut_WF_Snw;

for i=1:length(uo)
    T_WF_Snw(i,:)=1./V_Snw(i,:);
    n_ER_WF_Snw(i)=T_WF_Snw(i,1)/(i_Sw_shock(i)-1);
    ER_WF_Snw(i,1:i_Sw_shock(i))=0:n_ER_WF_Snw(i):T_WF_Snw(i,1);
end


for i=1:length(uo)
    for k=1:length(Snw)
        if k<i_Sw_shock(i)
            ER_WF_Snw(i,k)=ER_WF_Snw(i,k);
            T_WF_Snw(i,k)=ER_WF_Snw(i,k);
        else
            ER_WF_Snw(i,k)=Snw(k)-(fw(i,k)-1)/dfds(i,k);
        end
    end
end


ER_WF_Sw=ER_WF_Snw;
ER_WF_Sw=ER_WF_Snw*(1-Sor);
save Oil_Recovery_Noemalized.mat ER_WF_Snw 
save Oil_Recovery.mat ER_WF_Sw
save Water_Cut.mat Wcut_WF_Snw % fw actually
%% Test and Result
%  plot the relative permeability curve

figure(1) % Relative Permeability Curve
plot(Sw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),Krw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),'b',...
    Sw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),Kro(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Relative Permeability','fontsize',16)
legend({'Krw','Kro'},'fontsize',12)

saveas(figure(1),'Relative Permeability Curve.tif')

figure(2) % Fractional Flow Curve and Tangent Line, Water Sturation
plot(Sw,fw(1,:),'b',Sw,fw(2,:),'r',Sw,fw(3,:),'g','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
hold on

for i=1:length(uo)
    S_shock(i,:)=[Siw,Sw_shock(i)];
    F_shock(i,:)=[0,fw_shock(i)];
    S_shock_end(i,:)=[Siw,ER_BT_Sw(i)];
    F_shock_end=[0,1];
end

plot(S_shock(1,:),F_shock(1,:),'b-',S_shock(2,:),F_shock(2,:),'r-',S_shock(3,:),F_shock(3,:),'g-',...
    S_shock_end(1,:),F_shock_end,'b--',S_shock_end(2,:),F_shock_end,'r--',S_shock_end(3,:),F_shock_end,'g--','linewidth',1)
hold off

legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(2),'Fractional Flow and shock.tif')

figure(3) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(Snw,fw(1,:),'b',Snw,fw(2,:),'r',Snw,fw(3,:),'g','linewidth',2)
axis([0 1 0 1])
xlabel('Snw, Normalized Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
% legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)
hold on

for i=1:length(uo)
    Sn_shock(i,:)=[0,Snw_shock(i)];
    F_shock(i,:)=[0,fw_shock(i)];
    Sn_shock_end(i,:)=[0,ER_BT_Snw(i)];
    F_shock_end=[0,1];
end
%Sn_shock=[0,Snw_shock];
%F_shock=[0,fw_shock];
%Sn_shock_end=[0,ER_BT_Snw];
%F_shock_end=[0,1];
plot(Sn_shock(1,:),F_shock(1,:),'b-',Sn_shock(2,:),F_shock(2,:),'r-',Sn_shock(3,:),F_shock(3,:),'g-',...
    Sn_shock_end(1,:),F_shock_end,'b--',Sn_shock_end(2,:),F_shock_end,'r--',Sn_shock_end(3,:),F_shock_end,'g--','linewidth',1)
hold off

legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(3),'Fractional Flow and shock, Normalized.tif')

figure(4) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
loglog(T_WF_Snw(1,:),ER_WF_Snw(1,:),'b',T_WF_Snw(2,:),ER_WF_Snw(2,:),'r',T_WF_Snw(3,:),ER_WF_Snw(3,:),'g','linewidth',2)
axis([0 5000 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(4),'Normalized Recovery Factor, mobility ratio.tif')

figure(5) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw(1,:),ER_WF_Snw(1,:),'b',T_WF_Snw(2,:),ER_WF_Snw(2,:),'r',T_WF_Snw(3,:),ER_WF_Snw(3,:),'g','linewidth',2)
axis([0 50 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(5),'Normalized Recovery Factor, mobility ratio, first 50 PV.tif')

figure(6) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw(1,:),Wcut_WF_Snw(1,:),'b',T_WF_Snw(2,:),Wcut_WF_Snw(2,:),'r',T_WF_Snw(3,:),Wcut_WF_Snw(3,:),'g','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Water Cut','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(6),'Normalized Water Cut, mobility ratio, first 5 PV.tif')

figure(7) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw(1,:),Ocut_WF_Snw(1,:),'b',T_WF_Snw(2,:),Ocut_WF_Snw(2,:),'r',T_WF_Snw(3,:),Ocut_WF_Snw(3,:),'g','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Oil Cut','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(7),'Oil Cut, mobility ratio.tif')

figure(8) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw(1,:),ER_WF_Snw(1,:),'b',T_WF_Snw(2,:),ER_WF_Snw(2,:),'r',T_WF_Snw(3,:),ER_WF_Snw(3,:),'g','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor, Normalized','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(8),'Normalized Recovery Factor, mobility ratio, 5 PV.tif')

figure(9) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
plot(T_WF_Snw(1,:),ER_WF_Sw(1,:),'b',T_WF_Snw(2,:),ER_WF_Sw(2,:),'r',T_WF_Snw(3,:),ER_WF_Sw(3,:),'g','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(9),'Recovery Factor, mobility ratio.tif')

figure(10) % Water Saturation Profile
t=0.3;
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Snw,'b',X_Snw(2,:),Snw,'r',X_Snw(3,:),Snw,'g','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, X_D','fontsize',16)
ylabel('Normalized Water Saturation','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(10),'Normalized water saturation profile, mobility ratio.tif')

figure(11) % Water Saturation Profile
t=0.3;
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Sw,'b',X_Snw(2,:),Sw,'r',X_Snw(3,:),Sw,'g','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, X_D','fontsize',16)
ylabel('Water Saturation','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(11),'Water saturation profile, mobility ratio.tif')

figure(12) % Fractional Flow Curve and Tangent Line, Normalized Water Sturation
loglog(T_WF_Snw(1,:),ER_WF_Sw(1,:),'b',T_WF_Snw(2,:),ER_WF_Sw(2,:),'r',T_WF_Snw(3,:),ER_WF_Sw(3,:),'g','linewidth',2)
axis([0 5000 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))],['M^0=',num2str(M(3))]},'location','best','fontsize',12)

saveas(figure(12),'Recovery factor, 5000 PV, mobility ratio.tif')
