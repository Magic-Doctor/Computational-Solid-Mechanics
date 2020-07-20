% -----------------------------*杆单元*----------------------------------
% -----------------------------*静力计算*--------------------------------
clc;
clear;
AE = [200*10^9;10^-4;10^3]; %单元特征数（第一行存储弹性模量，第二行存储截面面积，第三行存储质量密度）
NMN = [1 1 1 1 1 1 1 1 1]; %单元特征数类信息(按单元存储)
ME = [3 2 5 1 1 1 2 4 3;
      4 3 6 4 2 5 5 6 6]; %单元信息(每个单元的节点编号)
Coordinates = [0 0;1 0;1 1;0 1;1/3 1/2;2/3 1/2]; %坐标(一行为一个节点的x,y,z坐标，按节点顺序填写(平面问题不填z坐标))
Elength = eLength(ME,Coordinates);%每个单元的长度
NE = size(ME,2);%总单元数
NP = size(Coordinates,1);%总节点数 
NF = size(Coordinates,2);%节点自由度(根据输入坐标判断节点自由度)
NRR = zeros(NF,NP);%节点约束信息(0未约束，1有约束)
NRR(1,1) = 1;NRR(2,1) = 1;NRR(2,2) = 1;%节点的约束位置
P = zeros(NF,NP);%节点载荷信息
P(2,6) = -1000;%节点载荷(约束处的载荷在Cholesky分解计算时均当作0(主对角元置大数时行列不会划去依然存在))
IS = zeros(NE,NF*2);%节点全局编号信息
K = zeros(NP*NF);%总刚
stabK = zeros(NP*2);%总体几何矩阵(此程序只能计算平面稳定性问题）
L_d = zeros(2,NE);%局部单元节点位移
L_F = zeros(2,NE);%局部单元节点力 
G_F = zeros(NF*2,NE);%全局单元节点力
PP = zeros(NF,NP);%结构节点力(将不同单元的同一个节点的单元节点力相加，抵消相互作用力)
SG = zeros(1,NE);%单元应力（按单元存储）


%计算总刚
for i=1:NE
    E = AE(1,NMN(i));A = AE(2,NMN(i));
    %局部坐标系单刚
    L_Ke = E*A/Elength(i) * [1 -1;-1 1];
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    %全局坐标系单刚
    G_Ke = TK' * L_Ke * TK;
    %节点全局编号
    for m = 1:NF*2
        if m <= NF
            IS(i,m) = (ME(1,i) - 1) * NF + m;
        else
            IS(i,m) = (ME(2,i) - 1) * NF + m - NF;
        end
    end
    %集成总刚
    for m = 1:NF*2
        for n = 1:NF*2
            K(IS(i,m),IS(i,n)) = K(IS(i,m),IS(i,n)) + G_Ke(m,n);
        end
    end
    fprintf("第%d号单元局部单刚\n\n",i);disp(L_Ke);
    fprintf("第%d号单元转换阵\n\n",i);disp(TK);
    fprintf("第%d号单元全局单刚\n\n",i);disp(G_Ke);
end
fprintf("总刚度阵\n\n");disp(K);

%约束处理
for j = 1:NP
    for i = 1:NF
        if NRR(i,j) == 1
            IS_v = (j - 1) * NF + i;
            K(IS_v,IS_v) = 10^30;%主对角元置大数
        end
    end
end

%计算全局节点位移
G_d = Cholesky(K,P,NF);


%计算局部节点位移/局部节点力/单元应力/节点力
for i = 1:NE
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    if NF == 2
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4))];%i号单元全局节点位移（平面下）
    else 
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4)) G_d(IS(i,5)) G_d(IS(i,6))];%i号单元全局节点位移（空间下）
    end
    L_d(:,i) = TK * E_d';%i号局部单元节点位移
    E = AE(1,NMN(i));A = AE(2,NMN(i));
    L_Ke = E*A/Elength(i) * [1 -1;-1 1];
    L_F(:,i) = L_Ke * L_d(:,i);%i号局部单元节点力
    SG(i) = L_F(2,i) / A;%i号单元应力(第二个节点的局部节点力正负能反映拉正压负）
    G_F(:,i) = TK' * L_F(:,i);%i号单元全局节点力
    if NF == 2
        s1 = [G_F(1,i) G_F(2,i)];s2 = [G_F(3,i) G_F(4,i)];
    else
        s1 = [G_F(1,i) G_F(2,i) G_F(3,i)];s2 = [G_F(4,i) G_F(5,i) G_F(6,i)];
    end
    PP(:,ME(1,i)) = PP(:,ME(1,i)) + s1';%将i号单元第一个节点的单元节点力存入PP列阵
    PP(:,ME(2,i)) = PP(:,ME(2,i)) + s2';%将i号单元第二个节点的单元节点力存入PP列阵
end

%计算支座反力
FR = PP - P;

%输出支座反力
fprintf("约束反力FR\n\n");
for j = 1:NP
    for i = 1:NF
        if NRR(i,j) == 1 %根据节点约束信息来判断那些节点施加了约束
            if i == 1
                fprintf("节点%dx方向约束反力",j);
                disp(FR(i,j));
            elseif i == 2
                fprintf("节点%dy方向约束反力",j);
                disp(FR(i,j));
            else 
                fprintf("节点%dz方向约束反力",j);
                disp(FR(i,j));
            end
        end
    end
end

%输出局部单元节点力
fprintf("\n\n局部单元节点力L_F\n\n");
disp(L_F);%矩阵的第二行正负代表实际杆件拉压状态

%输出局部单元节点位移
fprintf("\n\n局部单元节点位移L_d\n\n");
disp(L_d);%矩阵的第二行正负代表实际位移的正负

% %输出单元节点力
% fprintf("\n\n局部单元节点力G_F\n\n");
% disp(G_F);

%输出节点位移
fprintf("\n\n全局单元节点位移G_d\n\n");
disp(G_d);

%输出单元内力
fprintf("\n\n单元内力\n\n");
disp(L_F(2,:));%单元内力其实就是局部单元节点力第二行的值

%输出单元应力
fprintf("\n\n单元应力SG\n\n");
disp(SG);


if NF == 3
    errordlg("此程序无法计算空间的稳定性问题");
else
% -----------------------------*稳定性计算*------------------------------
    assumF = 1;%默认假设的外力为1N
    F = L_F(2,:) * assumF/P(2,6);%假设外力下计算出的单元内力
    for i = 1:NE
        %转换阵
        stabTK = stabTransformation(ME(:,i),Coordinates,Elength(i));
        %局部坐标下几何矩阵
        L_stabKe = F(i)/Elength(i) * [0 0 0 0;0 1 0 -1;0 0 0 0;0 -1 0 1];
        %全局坐标下几何矩阵
        G_stabKe = stabTK' * L_stabKe * stabTK;
        %集成总体几何矩阵
        for m = 1:NF*2
            for n = 1:NF*2
                stabK(IS(i,m),IS(i,n)) = stabK(IS(i,m),IS(i,n)) + G_stabKe(m,n);
            end
        end
        fprintf("第%d号单元局部几何阵\n\n",i);disp(L_stabKe);
        fprintf("第%d号单元转换阵\n\n",i);disp(stabTK);
        fprintf("第%d号单元全局几何阵\n\n",i);disp(G_stabKe);
    end
    fprintf("总体几何矩阵\n\n");disp(stabK);

    %划行划列
    [dele_K,dele_stabK] = Delete(K,stabK,NF,NP,NRR);
    fprintf("划行划列后总刚\n\n");
    disp(dele_K);
    fprintf("划行划列后几何矩阵\n\n");
    disp(dele_stabK);
    lamda = eig(dele_K,dele_stabK);%临界比值
    fprintf("临界力参数为\n\n");
    disp(lamda);
    criticalF = lamda  * assumF; %临界力
    fprintf("临界力为\n\n");
    disp(criticalF);
end


% --------------------------*结构无阻尼自由振动计算*----------------------
M = zeros(NF*NP);%总体协调质量阵
%计算总体协调质量阵
for i=1:NE
    density = AE(3,NMN(i));
    A = AE(2,NMN(i));
    %局部坐标系质量阵
    L_Me = density*A*Elength(i)/6 * [2 1;1 2];
    TK = Transformation(ME(:,i),Coordinates,Elength(i));
    %全局坐标系质量阵
    G_Me = TK' * L_Me * TK;
    %集成总刚
    for m = 1:NF*2
        for n = 1:NF*2
            M(IS(i,m),IS(i,n)) = M(IS(i,m),IS(i,n)) + G_Me(m,n);
        end
    end
    fprintf("第%d号单元局部质量阵\n\n",i);disp(L_Me);
    fprintf("第%d号单元转换阵\n\n",i);disp(TK);
    fprintf("第%d号单元全局质量阵\n\n",i);disp(G_Me);
end
fprintf("总体协调质量阵\n\n");disp(M);

%划行划列
[dele_K,dele_M] = Delete(K,M,NF,NP,NRR);
fprintf("划行划列后质量矩阵\n\n");
disp(dele_M);
[vibMode,a] = eig(dele_K,dele_M);%a为自振频率的平方
fprintf("自振频率的平方\n");
disp(a);
w=sqrt(a);%自振频率
fprintf("自振频率\n");
disp(w);
% fprintf("振型\n");
% disp(vibMode);