function z=find_shock_Com(S,M)
global no nw
fw=1/(1+((1-S)^no/(S^nw))/M); % fw=1./(1+Kro.*uw./(Krw.*uo));
dfds=((fw^2)/M)*(((1-S)^no)/(S)^nw)*(no/(1-S)+nw/(S)); % derivative of water fractional dlow
deltafs=fw/S; % slope of fractional flow

z=dfds-deltafs;
end