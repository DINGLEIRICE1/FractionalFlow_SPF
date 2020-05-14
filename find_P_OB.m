function z=find_P_OB(S)
global M no nw Siw Sor deltafs_shock_p

Snwf=(S-Siw)/(1-Siw-Sor);
for i=1:length(Snwf)
   if Snwf(i)<=0
       Snwf(i)=eps;
   elseif Snwf(i)>=1
       Snwf(i)=1-eps;
   end
end
 
fw=1/(1+((1-Snwf)^no/((Snwf)^nw))/M); % fw=1./(1+Kro.*uw./(Krw.*uo));
fw_shock_p=(deltafs_shock_p)*S; % slope of fractional flow
z=fw-fw_shock_p;
end