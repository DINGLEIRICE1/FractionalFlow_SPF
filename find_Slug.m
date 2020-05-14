function [A,B]=find_Slug(V_p,Sw_vp,Fw_vp)

A=Fw_vp-V_p*Sw_vp;
B=Fw_vp/V_p-Sw_vp;

end