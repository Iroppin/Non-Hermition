t1=1;t2=1;gamma=4/3;L=100;N=2*L; H=zeros(2*L,2*L); 
    for n=1:L
        if n>1
            H(2*n-1,2*n-2)=t2;
        end
        H(2*n-1,2*n)=t1+gamma/2;
        H(2*n,2*n-1)=t1-gamma/2;
        if n<L
            H(2*n,2*n+1)=t2;
        end
    end
% [vector_NNarray,eigN1array]=eig(H,'vector');
% abs_eig=abs(eigN1array);
% con_number=abs_eig(1);
% 
% for k=1:length(abs_eig)
%     middle=abs_eig(k);
%     if middle<con_number
%         con_number=middle;
[vector_NNarray,              eigN1array]=eig(H,'vector');%加了vector返回的特征值就是向量不是对角矩阵
  %    ↑每列为一个特征向量        ↑特征值               
abs_eig=abs(eigN1array);
[sortresult,                                     number_array]=                     sort(abs_eig);
%↑排序之后的结果，不过没啥用，其实可以用~代替          ↑前者中的abseig在原来位置的序号       ↑(￢_￢)本代码中最重要的一步
zero_number=number_array(1:2);
zero_mode=vector_NNarray(:,zero_number(1));
figure;
plot(1:N,abs(zero_mode).^2);%注意1：N
figure;
    all_number=setdiff(1:N,zero_number);
                             %除去零模
select_number=randperm(length(all_number),8);
select_result=all_number(select_number);
hold on;
for i=1:8
    state=vector_NNarray(:,all_number(i));
    plot(1:N,abs(state).^2);
end










