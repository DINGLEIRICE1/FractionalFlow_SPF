function z=find_shock_p(S)
global M_p no nw Siw Sor

Snwf=(S-Siw)/(1-Siw-Sor);
for i=1:length(Snwf)
   if Snwf(i)<=0
       Snwf(i)=eps;
   elseif Snwf(i)>=1
       Snwf(i)=1-eps;
   end
end
 
fw_p=1/(1+((1-Snwf)^no/((Snwf)^nw))/M_p); % fw=1./(1+Kro.*uw./(Krw.*uo));
dfds_p=((fw_p^2)/M_p/(1-Siw-Sor))*(((1-Snwf)^no)/(Snwf)^nw)*(no/(1-Snwf)+nw/(Snwf)); % derivative of water fractional dlow
deltafs_p=fw_p/S; % slope of fractional flow

z=dfds_p-deltafs_p;
end