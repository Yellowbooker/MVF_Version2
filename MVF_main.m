%% Second Step: Modal-Based Vector Fitting (MVF) Technique
% Ref: Circuit Model Extraction of Parallel-Connected Dual-Passband
% Coupled-Resonator Filters (Dr.Zhao Ping)

% ======================= Clear Workspace =======================

clear
close all
warning off
% The warning is about the polt issue... ignore it


format short
% ======================= Loop Parameters =======================

X = 18; % Max number of iterations
err_ak_limit = -120;    % in dB
err_cp_limit = -120;	% in dB

% ======================= Filter Parameters =======================


N = 8; % Order of the filter
Nz = 4; % Number of the transmission zeros 
% If the extracted Tz is not correct,plz increase Nz and 
% create a new parasitic coupling to explain the additional Tz

% Example
CF = 12316*1e6; % Hz center mapping frequency
BW = 36*1e6; % Hz mapping bandwidth
FBW = BW/CF;

% % Set fitting frequency range in the lowpass range
w1 = -6;     
w2 = 6;



% ======================= Sampling data import =======================

[F_band,Origin_S11,Origin_S12,Origin_S22]=Reads2p();


Ns=length(F_band);% number of sampling data

Origin_dBS11(:,1)=20*log10(abs(Origin_S11(:,1)));
Origin_dBS12(:,1)=20*log10(abs(Origin_S12(:,1)));
Origin_dBS22(:,1)=20*log10(abs(Origin_S22(:,1)));
w_low=(F_band./CF-CF./F_band)/FBW;
s_low=1i*w_low;

[empt,i1]=min(abs(w_low-ones(Ns,1)*w1));
[empt,i2]=min(abs(w_low-ones(Ns,1)*w2));

% De_Embedding_2p is a function to fit S11 and S22 independently
% De_Embedding_2p_cp is a function to fit S11 and S22 with the same poles

[S11,S12,S22]=De_Embedding_2p(F_band,Origin_S11,Origin_S12,Origin_S22,N,CF,BW);
% [S11,S12,S22]=De_Embedding_2p_cp(F_band,Origin_S11,Origin_S12,Origin_S22,N,CF,BW);

% find the fitting range
Ns = i2 - i1 + 1;
S11 = S11(i1:i2);
S12 = S12(i1:i2);
S22 = S22(i1:i2);
s_low = s_low(i1:i2);
F_band = F_band(i1:i2);
Origin_dBS11 = Origin_dBS11(i1:i2);
Origin_dBS12 = Origin_dBS12(i1:i2);
Origin_dBS22 = Origin_dBS22(i1:i2);

Z0=1;
for index=1:Ns
    [Y11(index,1),Y12(index,1),Y22(index,1)]=...
        StoY(Z0,S11(index,1),S12(index,1),S22(index,1));
end
dBY11(:,1)=20*log10(abs(Y11(:,1)));
dBY12(:,1)=20*log10(abs(Y12(:,1)));

% ======================= Initial pole ak =======================

for index=1:N
    ak(index,1)=-1+2.1/(N-1)*(index-1);
    ak(index,1)=-0.01*abs(ak(index,1))+ak(index,1)*1i;
end

% ======================= Iterate Max X times =======================

A=zeros(N,N);
b=ones(N,1);A1=zeros(Ns,N+1);A2=zeros(Ns,N);A3=zeros(Ns,Nz+1);
M_Y11 = diag(Y11);
M_Y12 = diag(Y12);
M_Y22 = diag(Y22);


% Weight of the fitting
% W12 = diag(1./sqrt(abs(Y12)));
W12 = diag(1./(abs(Y12).^0.5));
% W12 = eye(N_s);

tmp_ak = ak;
cp_plot = zeros(N, X);
err_ctl = zeros(1, X);

for index=1:X
    
    A_den=ones(Ns,1);
    
    for index1 = 1:N
        A2(:,index1) = 1./(s_low - ones(Ns,1).*ak(index1));
        A_den(:,1) = A2(:,index1).*A_den(:,1);
    end
    
    for index1 = 1:Nz+1
        A3(:,index1) = (s_low.^(index1-1)).*A_den;
    end
 
    A1 = [A2,ones(Ns,1)];
    A = diag(ak);

% %================Use QR factorization======================================
% 
%  % Thin QR factorization
% [Q1, R1] = qr([A1 -1*M_Y11*A2], 0);
% [Q2, R2] = qr([W12*A3 -1*W12*M_Y12*A2], 0);
% [Q3, R3] = qr([A1 -1*M_Y22*A2], 0);
% 
% % Formulate the LS solution for c
% B1_2 = Q1'*Y11;    A1_2 = R1(end-N+1:end, end-N+1:end);    B1_3 = B1_2(end-N+1:end);
% B2_2 = Q2'*Y12;    A2_2 = R2(end-N+1:end, end-N+1:end);    B2_3 = B2_2(end-N+1:end);
% B3_2 = Q3'*Y22;    A3_2 = R3(end-N+1:end, end-N+1:end);    B3_3 = B3_2(end-N+1:end);
% left = [ A1_2; A2_2; A3_2 ];
% right = [ B1_3; B2_3; B3_3 ];
% 
% % Solve c
% cp = left\right;
% cp_plot(:,index)=10*log10(abs(cp));
% ak=eig(A-b*cp.');
% err_ctl(index) = sum(abs(ak - tmp_ak));
% if (10*log10(err_ctl(index)) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
%     break
% end
% tmp_ak = ak;
% %============================end===========================================

%============================ Original=====================================    
    left=[A1                zeros(Ns,Nz+1)   zeros(Ns,N+1)   -1*M_Y11*A2;
          zeros(Ns,N+1)    W12*A3            zeros(Ns,N+1)   -1*W12*M_Y12*A2;
          zeros(Ns,N+1)    zeros(Ns,Nz+1)   A1               -1*M_Y22*A2];
    right=[Y11;
           W12*Y12;
           Y22];
    Call=lsqminnorm(left,right); % This function can be replaced by \
    
    C=mat2cell(Call,[N+1,Nz+1,N+1,N],1);
    cp=C{4,1};
    cp_plot(:,index)=10*log10(abs(cp));
    ak=eig(A-b*cp.');
	err_ctl(index) = sum(abs(ak - tmp_ak));
	if (10*log10(err_ctl(index)) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
		break
	end
	tmp_ak = ak;
%==============================end=========================================
end

%==============================	Disp loop time	======
disp('Loop time:')
disp(index)

figure('name','MVF-error')
subplot(1,2,1);
plot([1:index],max(cp_plot(:,[1:index])),'b','Linewidth',1.5);
legend('cp', 'Location', 'NorthEast')
title('Cp error in Log')
grid on

subplot(1,2,2);
plot([1:index],10*log10(err_ctl([1:index])),'r','Linewidth',1.5);
legend('ak', 'Location', 'NorthEast')
title('ak error in Log')
grid on

C11=C{1,1};
C12=C{2,1};
C22=C{3,1};

% ======================= Calculate ExtrY =======================

ExtrY11 = A1*C11;
ExtrY12 = A3*C12;
ExtrY22 = A1*C22;


dBExtrY11(:,1)=20*log10(abs(ExtrY11(:,1))); dBExtrY12(:,1)=20*log10(abs(ExtrY12(:,1)));
dBExtrY22(:,1)=20*log10(abs(ExtrY22(:,1)));


figure('name','Y1')
plot(imag(s_low),dBY11,'y*',imag(s_low),dBY12,'g*',imag(s_low),dBExtrY11,'r',imag(s_low),dBExtrY12,'b');



% ======================= Find the R12 =======================
for index=1:N
    for index1=1:N-1
        if index1>=index
            flag=1;
        else
            flag=0;
        end
        trans2(1,index1)=ak(index,1)-ak(index1+flag,1);
    end
    R12(index,1)=polyval(flipud(C12),ak(index,1))/prod(trans2);
end       

% ======================= Find the Tzs =======================

Tz=roots(flipud(C12));
display(Tz);

% ======================= Y to TCM =======================

K11=C11(N+1,1);K22=C22(N+1,1);
M=zeros(N+2,N+2);
for index=1:N
    M(index+1,index+1)=1i*ak(index,1);
    if abs(C11(index,1))>=abs(C22(index,1))
        M(1,index+1)=sqrt(C11(index,1));
        M(N+2,index+1)=R12(index,1)/sqrt(C11(index,1));
    else if abs(C11(index,1))<abs(C22(index,1))
        M(N+2,index+1)=sqrt(C22(index,1));
        M(1,index+1)=R12(index,1)/sqrt(C22(index,1));
        end
    end
    M(index+1,1)=M(1,index+1);
    M(index+1,N+2)=M(N+2,index+1);
end
M(1,1)=K11/1i;
M(N+2,N+2)=K22/1i;
%% Third Step: Transform Matrix
% Ref: Microwave Filters for Communication Systems: Fundamentals, Design and Applications

for index1 = 1:length(Tz)
    for index2 = 1:length(Tz)-1
        if abs(imag(Tz(index2))) >= abs(imag(Tz(index2+1)))
            Tem = Tz(index2+1);
            Tz(index2+1) = Tz(index2);
            Tz(index2) = Tem;
        end
    end
end

% % ======================= folded coupling matrix =======================

foldedCM=to_foldedCM(N,M);

% % % ======================= wheel coupling matrix =======================
% 
% wheelCM=to_wheelCM(N,foldedCM);
% 
% % % ======================= CT coupling matrix =======================
% 
% [CTCM,flag]=to_CTCM(N,Tz.',wheelCM);
% 
% % % ======================= CQ coupling matrix =======================
%  
% CQCM=to_CQCM(flag,N,CTCM,length(Tz));
% % display(CQCM);
% 
% % % ======================= NextCT coupling matrix =======================
% [CQCM]=to_CTCM2(N,Tz(length(Tz)),CQCM);
% CQCM = normalizeCM(CQCM);
% display(real(CQCM));

% for index=1:N
%     Q(index,1)=-1/imag(CQCM(index+1,index+1))/FBW;
% end
% 
% figure('name','Q')
% bar(Q,'m');
% title('Unloaded Quality Factor')
% grid on

[ExtrS11,ExtrS12]=analyseCM(N,F_band,foldedCM,CF,FBW,Ns);
dBExtrS11(1,:)=20*log10(abs(ExtrS11(1,:)));
dBExtrS12(1,:)=20*log10(abs(ExtrS12(1,:)));
figure('name','S-Parameters');
plot(F_band,Origin_dBS11,'y',F_band,Origin_dBS12,'g',F_band,dBExtrS11,'r--',F_band,dBExtrS12,'b--','Linewidth',1.5);
legend('Measured S11', 'Measured S21','Fitted S11','Fitted S21', 'Location', 'NorthWest')
title('S Parameter in Mag')
xlim([F_band(1),F_band(Ns)]);
grid on

% ReExtrS11(1,:)=real(ExtrS11(1,:));
% ImExtrS11(1,:)=imag(ExtrS11(1,:));
% figure('name','S-Parameters');
% plot(F_band,real(S11),'y',F_band,imag(S11),'g',F_band,ReExtrS11(1,:),'r--',F_band,ImExtrS11(1,:),'b--','Linewidth',1.5);
% legend('Measured reS11', 'Measured imS11','Fitted reS11','Fitted imS21', 'Location', 'NorthWest')
% title('S Parameter in Mag')
% xlim([F_band(1),F_band(N_s)]);
% grid on