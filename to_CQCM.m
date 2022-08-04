function CQCM=to_CQCM(a,b,c,d)
flag=a;%--------标志�?
N=b;%---------order of filter
M=c;%---------N+2orderCTCM
nfz=d;
if flag==0
   disp('too many Tzs');
else
   if nfz<2
     disp('more than 2 plz!');%------Tz>2才可以提取CQ结构
     NumTz=0;
%      fprintf('存在 %d 个完全四元组（CQ）\n',NumTz);
   else
     NumTz=fix(nfz/2);%-------------取整
%      fprintf('存在 %d 个完全四元组（CQ）\n',NumTz);
   end
end
% ----------------------
for i=1:NumTz
            k=4*i;
            m=k+1;
            l=m+1;
            n=l;
            theta=atan(M(k,l)/(M(m,n)));
            R=eye(N+2);
            R(k,k)=cos(theta);
            R(m,m)=R(k,k);
            R(k,m)=-sin(theta);
            R(m,k)=-R(k,m);
            M=R*M*R.';
            k=k-1;
            m=k+1;
            l=m+1;
            n=l;
            theta=atan(M(k,l)/(M(m,n)));
            R=eye(N+2);
            R(k,k)=cos(theta);
            R(m,m)=R(k,k);
            R(k,m)=-sin(theta);
            R(m,k)=-R(k,m);
            M=R*M*R.';
end
CQCM=M;
