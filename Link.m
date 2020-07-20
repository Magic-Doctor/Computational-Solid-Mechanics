% -----------------------------*�˵�Ԫ*----------------------------------
% -----------------------------*��������*--------------------------------
clc;
clear;
AE = [200*10^9;10^-4;10^3]; %��Ԫ����������һ�д洢����ģ�����ڶ��д洢��������������д洢�����ܶȣ�
NMN = [1 1 1 1 1 1 1 1 1]; %��Ԫ����������Ϣ(����Ԫ�洢)
ME = [3 2 5 1 1 1 2 4 3;
      4 3 6 4 2 5 5 6 6]; %��Ԫ��Ϣ(ÿ����Ԫ�Ľڵ���)
Coordinates = [0 0;1 0;1 1;0 1;1/3 1/2;2/3 1/2]; %����(һ��Ϊһ���ڵ��x,y,z���꣬���ڵ�˳����д(ƽ�����ⲻ��z����))
Elength = eLength(ME,Coordinates);%ÿ����Ԫ�ĳ���
NE = size(ME,2);%�ܵ�Ԫ��
NP = size(Coordinates,1);%�ܽڵ��� 
NF = size(Coordinates,2);%�ڵ����ɶ�(�������������жϽڵ����ɶ�)
NRR = zeros(NF,NP);%�ڵ�Լ����Ϣ(0δԼ����1��Լ��)
NRR(1,1) = 1;NRR(2,1) = 1;NRR(2,2) = 1;%�ڵ��Լ��λ��
P = zeros(NF,NP);%�ڵ��غ���Ϣ
P(2,6) = -1000;%�ڵ��غ�(Լ�������غ���Cholesky�ֽ����ʱ������0(���Խ�Ԫ�ô���ʱ���в��Ữȥ��Ȼ����))
IS = zeros(NE,NF*2);%�ڵ�ȫ�ֱ����Ϣ
K = zeros(NP*NF);%�ܸ�
stabK = zeros(NP*2);%���弸�ξ���(�˳���ֻ�ܼ���ƽ���ȶ������⣩
L_d = zeros(2,NE);%�ֲ���Ԫ�ڵ�λ��
L_F = zeros(2,NE);%�ֲ���Ԫ�ڵ��� 
G_F = zeros(NF*2,NE);%ȫ�ֵ�Ԫ�ڵ���
PP = zeros(NF,NP);%�ṹ�ڵ���(����ͬ��Ԫ��ͬһ���ڵ�ĵ�Ԫ�ڵ�����ӣ������໥������)
SG = zeros(1,NE);%��ԪӦ��������Ԫ�洢��


%�����ܸ�
for i=1:NE
    E = AE(1,NMN(i));A = AE(2,NMN(i));
    %�ֲ�����ϵ����
    L_Ke = E*A/Elength(i) * [1 -1;-1 1];
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    %ȫ������ϵ����
    G_Ke = TK' * L_Ke * TK;
    %�ڵ�ȫ�ֱ��
    for m = 1:NF*2
        if m <= NF
            IS(i,m) = (ME(1,i) - 1) * NF + m;
        else
            IS(i,m) = (ME(2,i) - 1) * NF + m - NF;
        end
    end
    %�����ܸ�
    for m = 1:NF*2
        for n = 1:NF*2
            K(IS(i,m),IS(i,n)) = K(IS(i,m),IS(i,n)) + G_Ke(m,n);
        end
    end
    fprintf("��%d�ŵ�Ԫ�ֲ�����\n\n",i);disp(L_Ke);
    fprintf("��%d�ŵ�Ԫת����\n\n",i);disp(TK);
    fprintf("��%d�ŵ�Ԫȫ�ֵ���\n\n",i);disp(G_Ke);
end
fprintf("�ܸն���\n\n");disp(K);

%Լ������
for j = 1:NP
    for i = 1:NF
        if NRR(i,j) == 1
            IS_v = (j - 1) * NF + i;
            K(IS_v,IS_v) = 10^30;%���Խ�Ԫ�ô���
        end
    end
end

%����ȫ�ֽڵ�λ��
G_d = Cholesky(K,P,NF);


%����ֲ��ڵ�λ��/�ֲ��ڵ���/��ԪӦ��/�ڵ���
for i = 1:NE
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    if NF == 2
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4))];%i�ŵ�Ԫȫ�ֽڵ�λ�ƣ�ƽ���£�
    else 
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4)) G_d(IS(i,5)) G_d(IS(i,6))];%i�ŵ�Ԫȫ�ֽڵ�λ�ƣ��ռ��£�
    end
    L_d(:,i) = TK * E_d';%i�žֲ���Ԫ�ڵ�λ��
    E = AE(1,NMN(i));A = AE(2,NMN(i));
    L_Ke = E*A/Elength(i) * [1 -1;-1 1];
    L_F(:,i) = L_Ke * L_d(:,i);%i�žֲ���Ԫ�ڵ���
    SG(i) = L_F(2,i) / A;%i�ŵ�ԪӦ��(�ڶ����ڵ�ľֲ��ڵ��������ܷ�ӳ����ѹ����
    G_F(:,i) = TK' * L_F(:,i);%i�ŵ�Ԫȫ�ֽڵ���
    if NF == 2
        s1 = [G_F(1,i) G_F(2,i)];s2 = [G_F(3,i) G_F(4,i)];
    else
        s1 = [G_F(1,i) G_F(2,i) G_F(3,i)];s2 = [G_F(4,i) G_F(5,i) G_F(6,i)];
    end
    PP(:,ME(1,i)) = PP(:,ME(1,i)) + s1';%��i�ŵ�Ԫ��һ���ڵ�ĵ�Ԫ�ڵ�������PP����
    PP(:,ME(2,i)) = PP(:,ME(2,i)) + s2';%��i�ŵ�Ԫ�ڶ����ڵ�ĵ�Ԫ�ڵ�������PP����
end

%����֧������
FR = PP - P;

%���֧������
fprintf("Լ������FR\n\n");
for j = 1:NP
    for i = 1:NF
        if NRR(i,j) == 1 %���ݽڵ�Լ����Ϣ���ж���Щ�ڵ�ʩ����Լ��
            if i == 1
                fprintf("�ڵ�%dx����Լ������",j);
                disp(FR(i,j));
            elseif i == 2
                fprintf("�ڵ�%dy����Լ������",j);
                disp(FR(i,j));
            else 
                fprintf("�ڵ�%dz����Լ������",j);
                disp(FR(i,j));
            end
        end
    end
end

%����ֲ���Ԫ�ڵ���
fprintf("\n\n�ֲ���Ԫ�ڵ���L_F\n\n");
disp(L_F);%����ĵڶ�����������ʵ�ʸ˼���ѹ״̬

%����ֲ���Ԫ�ڵ�λ��
fprintf("\n\n�ֲ���Ԫ�ڵ�λ��L_d\n\n");
disp(L_d);%����ĵڶ�����������ʵ��λ�Ƶ�����

% %�����Ԫ�ڵ���
% fprintf("\n\n�ֲ���Ԫ�ڵ���G_F\n\n");
% disp(G_F);

%����ڵ�λ��
fprintf("\n\nȫ�ֵ�Ԫ�ڵ�λ��G_d\n\n");
disp(G_d);

%�����Ԫ����
fprintf("\n\n��Ԫ����\n\n");
disp(L_F(2,:));%��Ԫ������ʵ���Ǿֲ���Ԫ�ڵ����ڶ��е�ֵ

%�����ԪӦ��
fprintf("\n\n��ԪӦ��SG\n\n");
disp(SG);


if NF == 3
    errordlg("�˳����޷�����ռ���ȶ�������");
else
% -----------------------------*�ȶ��Լ���*------------------------------
    assumF = 1;%Ĭ�ϼ��������Ϊ1N
    F = L_F(2,:) * assumF/P(2,6);%���������¼�����ĵ�Ԫ����
    for i = 1:NE
        %ת����
        stabTK = stabTransformation(ME(:,i),Coordinates,Elength(i));
        %�ֲ������¼��ξ���
        L_stabKe = F(i)/Elength(i) * [0 0 0 0;0 1 0 -1;0 0 0 0;0 -1 0 1];
        %ȫ�������¼��ξ���
        G_stabKe = stabTK' * L_stabKe * stabTK;
        %�������弸�ξ���
        for m = 1:NF*2
            for n = 1:NF*2
                stabK(IS(i,m),IS(i,n)) = stabK(IS(i,m),IS(i,n)) + G_stabKe(m,n);
            end
        end
        fprintf("��%d�ŵ�Ԫ�ֲ�������\n\n",i);disp(L_stabKe);
        fprintf("��%d�ŵ�Ԫת����\n\n",i);disp(stabTK);
        fprintf("��%d�ŵ�Ԫȫ�ּ�����\n\n",i);disp(G_stabKe);
    end
    fprintf("���弸�ξ���\n\n");disp(stabK);

    %���л���
    [dele_K,dele_stabK] = Delete(K,stabK,NF,NP,NRR);
    fprintf("���л��к��ܸ�\n\n");
    disp(dele_K);
    fprintf("���л��к󼸺ξ���\n\n");
    disp(dele_stabK);
    lamda = eig(dele_K,dele_stabK);%�ٽ��ֵ
    fprintf("�ٽ�������Ϊ\n\n");
    disp(lamda);
    criticalF = lamda  * assumF; %�ٽ���
    fprintf("�ٽ���Ϊ\n\n");
    disp(criticalF);
end


% --------------------------*�ṹ�����������񶯼���*----------------------
M = zeros(NF*NP);%����Э��������
%��������Э��������
for i=1:NE
    density = AE(3,NMN(i));
    A = AE(2,NMN(i));
    %�ֲ�����ϵ������
    L_Me = density*A*Elength(i)/6 * [2 1;1 2];
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    %ȫ������ϵ������
    G_Me = TK' * L_Me * TK;
    %�����ܸ�
    for m = 1:NF*2
        for n = 1:NF*2
            M(IS(i,m),IS(i,n)) = M(IS(i,m),IS(i,n)) + G_Me(m,n);
        end
    end
    fprintf("��%d�ŵ�Ԫ�ֲ�������\n\n",i);disp(L_Me);
    fprintf("��%d�ŵ�Ԫת����\n\n",i);disp(TK);
    fprintf("��%d�ŵ�Ԫȫ��������\n\n",i);disp(G_Me);
end
fprintf("����Э��������\n\n");disp(M);

%���л���
[dele_K,dele_M] = Delete(K,M,NF,NP,NRR);
fprintf("���л��к���������\n\n");
disp(dele_M);
[vibMode,a] = eig(dele_K,dele_M);%aΪ����Ƶ�ʵ�ƽ��
fprintf("����Ƶ�ʵ�ƽ��\n");
disp(a);
w=sqrt(a);%����Ƶ��
fprintf("����Ƶ��\n");
disp(w);
% fprintf("����\n");
% disp(vibMode);