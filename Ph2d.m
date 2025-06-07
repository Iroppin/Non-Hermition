L=40;t2=1;gamma=4/3;t1_range=-3:0.01:3;
points=length(t1_range);
molE=cell(points,1);
disturb=0.8;
for i=1:points
    t1=t1_range(i);
    H=zeros(2*L,2*L);
for n=1:L
    if n>1
        H(2*n-1,2*n-2)=t2;
        H(2*n-1,2*n)=t1+gamma/2;
        H(2*n,2*n-1)=t1-gamma/2;
    end
        H(1,2)=(t1-disturb)+gamma/2;
        H(2,1)=(t1-disturb)-gamma/2;
        if n<L
            H(2*n,2*n+1)=t2;
        end
end
molE{i}=eig(H);
end
% 就改一改H的写法，没区别

figure;
hold on;
for k=1:points
    t1=t1_range(k);
    E_p=molE{k};
    scatter(t1*ones(size(E_p)),abs(E_p),1,'k');
    zero_modes=E_p(abs(E_p)<0.1);
    scatter(t1*ones(size(zero_modes)),abs(zero_modes),10,'r');
end
hold off;
