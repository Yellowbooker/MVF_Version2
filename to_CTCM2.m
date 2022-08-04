function [CTCM]=to_CTCM2(a,b,c)
N=a;%------order of filter
Tz=b./1i;%---------location of Tz
nfz=length(b);
M=c;%----------N+2orderWheelCM
%--------------------------------CT,afterWheel

if (N/(nfz*2))==1
    NumTz=nfz-1;
    flag=1;
elseif (N/(nfz*2))>1
    NumTz=nfz;
    flag=1;
else
    disp('ERROR');
    NumTz=0;
    flag=0;
end

theta=atan(M(N,N+1)/(M(N+1,N+1)+Tz(1,1)));
R=eye(N+2);
R(N,N)=cos(theta);
R(N+1,N+1)=R(N,N);
R(N,N+1)=-sin(theta);
R(N+1,N)=-R(N,N+1);
M=R*M*R.';
CTCM=M;
