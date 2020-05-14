function z=find_shock_SP_Ads(S)
global M_SP no_SF nw_SF Siw_SF Sor_SF Retar_Norm

Snwf=(S-Siw_SF)/(1-Siw_SF-Sor_SF);
for i=1:length(Snwf)
   if Snwf(i)<=0
       Snwf(i)=eps;
   elseif Snwf(i)>=1
       Snwf(i)=1-eps;
   end
end
 
fw_SF=1/(1+((1-Snwf)^no_SF/((Snwf)^nw_SF))/M_SP); % fw=1./(1+Kro.*uw./(Krw.*uo));
dfds_SF=((fw_SF^2)/M_SP/(1-Siw_SF-Sor_SF))*(((1-Snwf)^no_SF)/(Snwf)^nw_SF)*(no_SF/(1-Snwf)+nw_SF/(Snwf)); % derivative of water fractional dlow
deltafs_SF=fw_SF/(S+Retar_Norm); % slope of fractional flow

z=dfds_SF-deltafs_SF;
end