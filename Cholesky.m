%-----------------------------*�ڵ�λ��*----------------------------------
function [G_d,L] = Cholesky(K,P,NF)%PΪNF��NP����
n = length(K);
L = zeros(n);
G_d = zeros(1,n);
P1 = zeros(1,n);

%��NF��NP��P����ת��Ϊ1��NF*NP����
for j = 1:n/NF %�ڵ���
    for i = 1:NF
        IS_v = (j - 1) * NF + i;
        P1(IS_v) = P(i,j);
    end
end

% ����˹���ֽ�
L(1,1) = sqrt(K(1,1));
for i = 2:n
    L(1,i) = K(1,i) / L(1,1);
    L(i,1) = L(1,i);
end

for i = 2:n
    for j = i:n
        a = 0;
        for k = 1:i-1
            a = a + L(i,k) * L(k,j);
        end
        if j == i 
            L(i,j) = sqrt(K(i,j) - a);
        else
            L(i,j) = (K(i,j) - a) / L(i,i);
            L(j,i) = L(i,j);
        end
    end
end

% ��Lx=P1
G_d(1) = P1(1,1) / L(1,1);
for i = 2:n
    a = 0;
    for j = 1:i-1
        a = a + L(i,j) * G_d(j);
    end
    G_d(i) = (P1(i) - a) / L(i,i);
end

%���L'x=x
G_d(n) = G_d(n) / L(n,n);
for i = n-1:-1:1
    a = 0;
    for j = n:-1:i+1
        a = a + L(i,j) * G_d(j);
    end
    G_d(i) = (G_d(i) - a) / L(i,i);
end