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
function [New_S11,New_S12,New_S22]=De_Embedding_2p(F_band,S11,S12,S22,N,CF,BW)
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
Sxx{1} = S11;
Sxx{2} = S22;

for INDEX = 1:length(Sxx)
    
    M_Sxx{INDEX} = diag(Sxx{INDEX});
    
    % ======================= Initial Common Pole ak [-3,2.3]

    for index = 1:N+M
        ak(index,1) = -3 + 5.3/(N + M - 1)*(index - 1);
        ak(index,1) = -0.01*abs(ak(index,1)) + ak(index,1)*1i;
    end
    tmp_ak = ak;
    cp_plot = zeros(N+M, X);
    err_ctl = zeros(1, X);

    % ======================= Fit S11 & S22

    for index = 1:X
        
        for index1 = 1:N+M
            A2(:,index1) = 1./(s_low - ones(N_s,1).*ak(index1));
        end
        
        A1 = [A2,ones(N_s,1)];
        A = diag(ak,0);
        left = [A1   -1*M_Sxx{INDEX}*A2];
        right = Sxx{INDEX};
    %     Call = left\right;
        C_all = lsqminnorm(left,right);
        C = mat2cell(C_all,[N+M,1,N+M],1);
        cp = C{3,1};
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
    pole{INDEX} = ak;
    residue{INDEX} = C{1,1};
    d{INDEX} = C{2,1};
    A11 = diag(pole{INDEX},0);
    zero{INDEX} = eig(A11-b*residue{INDEX}.'/d{INDEX});
    % ======================= Caculate Phase Factor
    index2 = 0;
    index3 = 0;
    for index=1:N+M
            if abs(zero{INDEX}(index,1))>=Line
                    index2 = index2+1;
                    alfa_zk{INDEX}(index2,1) = zero{INDEX}(index,1);
            end
            if abs(pole{INDEX}(index,1))>=Line
                    index3 = index3+1;
                    alfa_ak{INDEX}(index3,1) = pole{INDEX}(index,1);
            end
    end
    if index2 == 0 && index3 == 0
        for index = 1:N_s
        Extralfa{INDEX}(index,1) =1;
        end
    else
    TestNum = poly(alfa_zk{INDEX});
    TestDen = poly(alfa_ak{INDEX});
    for index = 1:N_s
        Extralfa{INDEX}(index,1) = d{INDEX}*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
        ExtrSxx{INDEX}(index,1) = d{INDEX}*polyval(poly(zero{INDEX}),s_low(index,1))/polyval(poly(pole{INDEX}),s_low(index,1));
    end
    end
    unwrap_Palfa{INDEX} = unwrap(angle(Extralfa{INDEX}));
    
end

figure('name','P-Z of S11')
plot(real(zero{1}),imag(zero{1}),'ro',real(pole{1}),imag(pole{1}),'bx','Linewidth',1.5);
legend('Zero S11', 'Pole S11', 'NorthWest')
figure('name','P-Z of S22')
plot(real(zero{2}),imag(zero{2}),'ro',real(pole{2}),imag(pole{2}),'bx','Linewidth',1.5);
legend('Zero S22', 'Pole S22', 'NorthWest')
figure('name','S11&ES11')
plot(w_low,abs(S11),'y',w_low,abs(ExtrSxx{1}),'b--','Linewidth',1.5);
figure('name','S22&ES22')
plot(w_low,abs(S22),'y',w_low,abs(ExtrSxx{2}),'b--','Linewidth',1.5);
% ======================= Unwrap Phase In Matlab Three Ports

New_S11 = S11./abs(Extralfa{1}).*exp(-1i*unwrap_Palfa{1});
New_S22 = S22./abs(Extralfa{2}).*exp(-1i*unwrap_Palfa{2});

New_S12 = S12./sqrt(abs(Extralfa{1}))./sqrt(abs(Extralfa{2})).*exp(-1i*(unwrap_Palfa{1}*0.5+unwrap_Palfa{2}*0.5));

figure('name','After')
plot(w_low,angle(New_S11),'r',w_low,angle(New_S22),'b','LineWidth',1.5);
title('S11&S22 in Phase(after de embedding)')
legend('S11', 'S22', 'NorthWest')
grid on

% ======================= Polt Unwraped Phase

figure('name','Phase loading')
subplot(1,2,1);
yyaxis left
plot(w_low,abs(Extralfa{1}),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa{1},'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S11')
grid on
subplot(1,2,2);
yyaxis left
plot(w_low,abs(Extralfa{2}),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa{2},'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S22')
grid on
