function [DeleMatrix1,DeleMatrix2] = Delete(Matrix1,Matrix2,NF,NP,NRR)
DeleMatrix1 = Matrix1;
DeleMatrix2 = Matrix2;
count = [];%�洢Լ�����ţ��Ա�ͬʱɾ�����ж���
for j = 1:NP %����Լ����
    for i = 1:NF 
        if NRR(i,j) == 1
            count = [count,(j-1)*NF+i];%�𲽵İ�Լ�����Ŵ洢��count_K����
        end
    end
end
DeleMatrix1(count,:) = [];%����Matrix1����
DeleMatrix1(:,count) = [];%����Matrix1����
DeleMatrix2(count,:) = [];%����Matrix2����
DeleMatrix2(:,count) = [];%����Matrix2����