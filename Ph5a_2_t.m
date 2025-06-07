function t3_nonzero
    % 主函数：计算t3≠0时的拓扑不变量和广义布里渊区
    tic;
    
    %参数
    t1=1.1; tx=1.1;
    t2=1;
    t3=0.2;  
    gamma=4/3;
    L=100;      
    delta_1=1e-8; % 微分步长
    delta_2=1e-3; % 积分步长。这两个步长都应该设置的比较短比较好，在t3为零的情况下，相变点在1.2+左右，且本应该竖直的直线有一个肉眼不易直接看出的斜率。就能看出来精读不够，但是此处精读的修改对运行时间的影响极大，在积分与微分的精读都上升一个数量级后，我跑了接近一个上午都没跑完。这里的参数目前是为了方便验证修改的参数，五秒之内就能跑完。如果要利用该代码复现文献的结果的话，需要适当提高一些精度。
    tol=0.01;     % β根筛选容差
    
    % 计算Cβ
    [beta_list,E]=calculate_beta_trajectory(t1,t2,t3,gamma,L,tol);
    
    
    % 计算绕数
    t1_range=linspace(-3,3,201); 
    W_values=zeros(size(t1_range));
    
    parfor i=1:length(t1_range)
        W_values(i)=calculate_winding_number(t1_range(i),t2,t3,gamma,delta_1,delta_2);
    end
    
    plot_winding_number(t1_range,W_values); % 绘制图5(a)
    
    toc;%计时用的，可以删掉
end

function [beta_list,E]=calculate_beta_trajectory(t1,t2,t3,gamma,L,tol)
    % 构造非厄米SSH模型的实空间哈密顿量并求解β根
    H=build_real_space_Hamiltonian(t1,t2,t3,gamma,L);
    E=eig(H); % 获取所有本征值
    
    % 求解每个能量E对应的β根
    beta_list=[];
    a=t1+gamma/2;
    b=t1-gamma/2;
    
    for e=E.'
        % 特征方程系数（四次方程）
        c4=t3*t2;
        c3=a*t2+t3*b;
        c2=t2^2+a*b+t3^2-e^2;
        c1=t2*b+a*t3;
        c0=t3*t2;
        
        % 求解四次方程根
        roots_beta=roots([c4,c3,c2,c1,c0]);
        
        % 筛选满足|βj|=|βk|的根对
        for j=1:4
            for k=j+1:4
                if abs(abs(roots_beta(j))-abs(roots_beta(k)))<tol
                    beta_list=[beta_list;roots_beta([j,k])];
                end
            end
        end
    end
end

function H=build_real_space_Hamiltonian(t1,t2,t3,gamma,L)
    % 构造实空间哈密顿量（包含t3项）
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
end

function W=calculate_winding_number(t1,t2,t3,gamma,delta_1,delta_2)
    % 计算非布洛赫绕数
    k=-pi:delta_2:pi-delta_2;
    k_delta=k+delta_1;
    W_integral=0;
    %下面也许需要Parallel Computing Toolbox
    for idx=1:length(k)
        beta=solve_beta_for_k(k(idx),t1,t2,t3,gamma);
        beta_delta=solve_beta_for_k(k_delta(idx),t1,t2,t3,gamma);
        
        % 计算H[0,1]项的对数导数，计算对数差值避免循环
        H01=(t1-gamma/2)+t2*beta+t3/beta;
        H01_delta=(t1-gamma/2)+t2*beta_delta+t3/beta_delta;
        
        d_log=(log(H01_delta)-log(H01))/delta_1;
        W_integral=W_integral+d_log*delta_2;
    end
    
    W=imag(W_integral)/(2*pi);
end

function beta=solve_beta_for_k(k,t1,t2,t3,gamma)
    % 数值求解广义布里渊区β
    % 此处简化为使用近似解析解，更严谨来说应该使用迭代。但是数值复现来说完全足够了
    a=t1+gamma/2;
    b=t1-gamma/2;
    
    % 近似处理
    beta=sqrt(abs((t1-gamma/2)/(t1+gamma/2)))*exp(1i*k); 
end


function plot_winding_number(t1_range,W_values)
    figure;
    

    % 自适应平滑窗口函数
    %也许是精度的原因，也许需要用文献中矩阵处理方法才行（用文献的方法我总是复现不了，遂如前所述），运行结果需要取整/降燥才能比较好的复现出论文的结果，这里简单使用平滑。其实在筛选beta的阈值那里也可以进行适当调整，但是阈值在0.1左右的时候才能较好的对beta轨迹进行复现，阈值再低的话很多的beta点都会取不到。实际上用Signal
    %Processing
    %Toolbox工具箱处理效果比较方便，这里仅仅简单进行平滑处理，但理论上高精度下不处理也没有什么问题
    window=0.5*(tanh(100*(t1_range))-tanh(100*(t1_range-1.1)));
    W_corrected=round(real(W_values)).*(1-window)+1.*window;
    
    %绘图
    plot(t1_range,W_corrected,'LineWidth',2);
    grid on;
    ylim([-0.5 1.5]);
    
  
end