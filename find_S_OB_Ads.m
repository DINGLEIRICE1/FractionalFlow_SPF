function z=find_S_OB_Ads(S)
global M no nw Siw Sor deltafs_shock_SF Retar_Norm

Snwf=(S-Siw)/(1-Siw-Sor);
for i=1:length(Snwf)
   if Snwf(i)<=0
       Snwf(i)=eps;
   elseif Snwf(i)>=1
       Snwf(i)=1-eps;
   end
end
 
fw=1/(1+((1-Snwf)^no/((Snwf)^nw))/M); % fw=1./(1+Kro.*uw./(Krw.*uo));
fw_shock_SF=(deltafs_shock_SF)*(S+Retar_Norm); % slope of fractional flow
z=fw-fw_shock_SF;
end