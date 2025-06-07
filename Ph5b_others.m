
t1=1.1;t2=1;gamma=4/3;t3=0.2;theta=linspace(0,2*pi,1000);%参数设置
r=zeros(size(theta));
options=optimset('Display','off','TolX',1e-6,'TolFun',1e-6);%使用工具箱

%数值求解每个角度对应的r
for i=1:length(theta)
    th=theta(i);
    %定义方程：|A|^2=|B|^2，其中A和B为特征方程的两个因子
    A=@(r_val)t1-gamma/2+t2*r_val*exp(1i*th)+t3/(r_val*exp(1i*th));
    B=@(r_val)t1+gamma/2+t2/(r_val*exp(1i*th))+t3*r_val*exp(1i*th);
    fun=@(r_val)abs(A(r_val))^2-abs(B(r_val))^2;
    %初始猜测：基于t3=0的解析解，并调整以考虑t3的影响
    r_initial=1;
    try
        r_sol=fzero(fun,r_initial,options);
        r(i)=r_sol;
    catch

    end
end
beta=r_valid.*exp(1i*theta_valid);
%绘制Cβ
figure;
scatter(real(beta),imag(beta),10,'filled');
axis equal;
grid on;