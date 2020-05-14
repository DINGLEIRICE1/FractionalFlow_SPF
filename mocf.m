%     FIND TANGENT FOR MOC.M
function z=asg12f(S)
global M n1 n2
i=1;
f=1 /(1+((1-S)^n2./(S^n1+eps))./(M(i)+eps));
dfds=(f^2 /M).*((1-S+eps)^n2/(S+eps)^n1)...
  *(n2/(1-S+eps)+n1/(S+eps));
z=f/(S+eps)-dfds;
