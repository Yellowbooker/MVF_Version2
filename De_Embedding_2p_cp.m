%% Microwave Filter Coupled-Resonator Circuit Modal Extraction By MVF
%% Version01----03/09/2020
%% Version02----03/27/2020
%% @Useless_Crap
%% First Step: Phase De-Embedding By VF
% Ref: Phase De-Embedding of Narrowband CoupledResonator Networks by
% Vector Fitting (Dr.Zhao Ping)

% % ======================= Clear Workspace  =======================
% 
% clear all
% clc
% 
% % ======================= Loop Parameters =======================
function [New_S11,New_S12,New_S22]=De_Embedding_2p_cp(F_band,S11,S12,S22,N,CF,BW)
X = 18;
err_ak_limit = -60;  %in dB
err_cp_limit = -60;	% in dB

% ======================= Phase Corrected  =======================

% ======================= Define Filter Parameter

FBW = BW/CF;

% ======================= Define Loop Parameter

M = 3; % Additional order to fit the phase loading (3 is typical)
Line = 2; % Radius that distinguish the poles and zeros (2 is typical)

% ======================= Bandpass to Lowpass

N_s = length(F_band); % number of sampling data

w_low = (F_band./CF - CF./F_band)/FBW;
s_low = 1i*w_low;

A = zeros(N+M,N+M);
b = ones(N+M,1);
A1 = zeros(N_s,N+M+1);
A2 = zeros(N_s,N+M);

M_S11 = diag(S11);
M_S22 = diag(S22);

% ======================= Initial Common Pole ak [-3,2.3]

for index = 1:N+M
    ak(index,1) = -3 + 5.3/(N + M - 1)*(index - 1);
    ak(index,1) = -0.01*abs(ak(index,1)) + ak(index,1)*1i;
end
tmp_ak = ak;
cp_plot = zeros(N+M, X);
err_ctl = zeros(1, X);

% ======================= Fit S11

for index = 1:X

    for index1 = 1:N+M
        A2(:,index1) = 1./(s_low - ones(N_s,1).*ak(index1));
    end

    A1 = [A2,ones(N_s,1)];
    A = diag(ak,0);
    left = [A1                zeros(N_s,N+M+1)     -1*M_S11*A2;
            zeros(N_s,N+M+1)  A1                   -1*M_S22*A2];
    right = [S11;
             S22];
%     Call = left\right;
    C_all = lsqminnorm(left,right);
    C = mat2cell(C_all,[N+M,1,N+M,1,N+M],1);
    cp = C{5,1};
    cp_plot(:,index)=10*log10(abs(cp));
    ak = eig(A-b*cp.');
    err_ctl(index) = sum(abs(ak - tmp_ak));
    if 10*log10(err_ctl(index)) < err_ak_limit && max(cp_plot(:,index)) < err_cp_limit
        break
    end
    tmp_ak = ak;
end


disp('Loop time:')
disp(index)


pole = ak;
residue11 = C{1,1};
d11 = C{2,1};
A = diag(pole,0);
zero11 = eig(A-b*residue11.'/d11);
residue22 = C{3,1};
d22 = C{4,1};
zero22 = eig(A-b*residue22.'/d22);

% ======================= Caculate Phase Factor
index3 = 0;
index2 = 0;
index1 = 0;
for index=1:N+M
        if abs(zero11(index,1))>=Line
                index1 = index1+1;
                alfa_zk11(index1,1) = zero11(index,1);
        end
        if abs(zero22(index,1))>=Line
                index2 = index2+1;
                alfa_zk22(index2,1) = zero22(index,1);
        end
        if abs(pole(index,1))>=Line
                index3 = index3+1;
                alfa_ak(index3,1) = pole(index,1);
        end
end

if index1 == 0 && index2 == 0 &&  index2 == 0
    for index = 1:N_s
    Extralfa11(index,1) =1;
    Extralfa22(index,1) =1;
    end
else
    TestNum11 = poly(alfa_zk11);
    TestNum22 = poly(alfa_zk22);
    TestDen = poly(alfa_ak);
    for index = 1:N_s
        Extralfa11(index,1) = d11*polyval(TestNum11,s_low(index,1))/polyval(TestDen,s_low(index,1));
        ExtrS11(index,1) = d11*polyval(poly(zero11),s_low(index,1))/polyval(poly(pole),s_low(index,1));
        Extralfa22(index,1) = d22*polyval(TestNum22,s_low(index,1))/polyval(TestDen,s_low(index,1));
        ExtrS22(index,1) = d22*polyval(poly(zero22),s_low(index,1))/polyval(poly(pole),s_low(index,1));
    end
end

unwrap_Palfa11 = unwrap(angle(Extralfa11));
unwrap_Palfa22 = unwrap(angle(Extralfa22));   

figure('name','P-Z of S11')
plot(real(zero11),imag(zero11),'ro',real(pole),imag(pole),'bx','Linewidth',1.5);
legend('Zero S11', 'Pole S11', 'NorthWest')
figure('name','P-Z of S22')
plot(real(zero22),imag(zero22),'ro',real(pole),imag(pole),'bx','Linewidth',1.5);
legend('Zero S22', 'Pole S22', 'NorthWest')
% figure('name','S11&ES11')
% plot(w_low,abs(S11),'y',w_low,abs(ExtrS11),'b--','Linewidth',1.5);
% figure('name','S22&ES22')
% plot(w_low,abs(S22),'y',w_low,abs(ExtrS22),'b--','Linewidth',1.5);
% ======================= Unwrap Phase In Matlab Three Ports

New_S11 = S11./abs(Extralfa11).*exp(-1i*unwrap_Palfa11);
New_S22 = S22./abs(Extralfa22).*exp(-1i*unwrap_Palfa22);

New_S12 = S12./sqrt(abs(Extralfa11))./sqrt(abs(Extralfa22)).*exp(-1i*(unwrap_Palfa11*0.5+unwrap_Palfa22*0.5));

figure('name','After')
plot(w_low,angle(New_S11),'r',w_low,angle(New_S22),'b','LineWidth',1.5);
title('S11&S22 in Phase(after de embedding)')
legend('S11', 'S22', 'NorthWest')
grid on

% ======================= Polt Unwraped Phase

figure('name','Phase loading')
subplot(1,2,1);
yyaxis left
plot(w_low,abs(Extralfa11),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa11,'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S11')
grid on
subplot(1,2,2);
yyaxis left
plot(w_low,abs(Extralfa22),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa22,'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S22')
grid on
