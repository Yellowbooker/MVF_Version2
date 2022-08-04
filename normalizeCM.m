function M = normalizeCM(M)
% Ensure main coupling larger than zero

N = length(M);


for index = 1 : N-1
    if (M(index, index+1) < 0)
        IM = eye(N);
        IM(index+1, index+1) = -1;
        M = IM * M * IM.';
    end
end

for index1 = 1 : N
    for index2 = 1 : N
        if (abs(M(index1, index2)) < 0.00001)
            M(index1, index2) = 0;
        end
    end
end
end