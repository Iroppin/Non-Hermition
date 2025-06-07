t1=1.1;
t2=1;
gamma=4/3;
t3=0.0001;
L=100;
H=zeros(2*L,2*L);
for n=1:L 
H(2*n-1,2*n)=t1+gamma/2;
if n>1
H(2*n-1,2*(n-1))=t2;
end
if n<L
H(2*n-1,2*(n+1))=t3;
end
H(2*n,2*n-1)=t1-gamma/2;
if n>1
H(2*n,2*(n-1)-1)=t3;
end
if n<L
H(2*n,2*(n+1)-1)=t2;
end
end
E=eig(H);
beta_list = [];
tol=0.01; 
for i=1:length(E)
e=E(i);
a=t1+gamma/2;
b=t1-gamma/2;
%最后还是选择直接解筛根，一开始我就用的这个思路，但是筛选步骤略有不慎……
c4=t3*t2;
c3=a*t2+t3*b;
c2=t2^2+a*b+t3^2-e^2;
c1=t2*b+a.*t3;
c0=t3*t2;
roots_beta=roots([c4,c3,c2,c1,c0]);

for j=1:4
beta_j=roots_beta(j);
for k=j+1:4%这种写法可以极大的简化循环中的结构
beta_k=roots_beta(k);
if abs(abs(beta_j)-abs(beta_k))<tol
beta_list=[beta_list; beta_j;beta_k];
end
end
end
end
figure;
scatter(real(beta_list), imag(beta_list),'filled');
axis equal;
grid on;