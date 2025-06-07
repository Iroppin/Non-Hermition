L=35;t2=1;gamma=4/3;t1_range=-3:0.01:3;
points=length(t1_range);
molE=cell(points,1);

for i=1:points
    t1=t1_range(i);
    H=zeros(2*L,2*L);
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
molE{i}=eig(H);

end
% 这里eig（H)是个2L*1的矩阵，molE的每个格子放一个不同t1下的本征值集合
figure;
subplot(2,2,1);
hold on;
for k=1:points
    t1=t1_range(k);
    E_p=molE{k};
    scatter(t1*ones(size(E_p)),abs(E_p),1,'k');
    %由于每个t1对应着80个E，需要ones对齐。
    zero_modes=E_p(abs(E_p)<0.1);
    scatter(t1*ones(size(zero_modes)),abs(zero_modes),10,'r');
end
hold off;

subplot(2,2,2);
hold on;
for k=1:points
    t1=t1_range(k);
    E_p=molE{k};
    scatter(t1*ones(size(E_p)),real(E_p),1,'k');
    zero_modes=E_p(abs(E_p)<0.1);
    scatter(t1*ones(size(zero_modes)),real(zero_modes),10,'r');
end
hold off;

subplot(2,2,3);
hold on;
for k=1:points
    t1=t1_range(k);
    E_p=molE{k};
    scatter(t1*ones(size(E_p)),imag(E_p),1,'k');
    zero_modes=E_p(abs(E_p)<0.1);
    scatter(t1*ones(size(zero_modes)),imag(zero_modes),10,'r');
end
hold off;







