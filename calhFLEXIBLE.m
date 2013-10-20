function H = calhFLEXIBLE(f,g,s,P,Q,pre1,pre2,w,cr,dof)
%计算弹性体的六向频响函数矩阵，如机组质心，浮筏，浮筏和机组的连接点之间，以及基座连接点之间。
%P和Q分别代表频响矩阵之中，激励点和响应点的个数。对于机组质心，P,Q=1；对于浮筏，P，Q可以分别为上下层连接点的数目，
%pre1和pre2可以取0，或者P和Q；对于机组的连接点之间，P，Q为每个机组下面的连接点的数目；对于基座，P，Q即为下层连接点的数目

%g的存储格式是：第1列为x向，第2列为y向....第6列为rz模态
    %存储格式是：
    %第1行：第一个连接点的第一阶模态位移
    %第2行：第1个连接点的第2阶模态位移
    %...
    %第N行：第1个连接点的第N阶模态位移
    %第N+1行：第2个连接点的第1阶模态位移
    %第N+2行：第2个连接点的第2阶模态位移
    %...
    %第2N行：第2个连接点的第N阶模态位移
        
[M,N] = size(f);%确定子结构频率的个数
H = zeros(P,Q);

for cn = 1 : dof    
    cnn = cn;
    mda = g(:,2); %取出第cnn列数据，不同方向
    
    for cm = 1 : dof        
        cmm = cm;
        mdb = s(:,2);   %取出第cmm列数据，不同方向
        
        for j=1:P/dof;%激励点数目
          for k=1:Q/dof;%响应点数目
              for mn=1:M;%频率              
              wr = 2*pi*f(mn,2);%角频率
              a = mda(pre1*M + (j-1)*M+mn,1);%a是第j点的第M阶模态在cn方向的模态阵型
              b = mdb(pre2*M + (k-1)*M+mn,1);
              
              H(dof*(j-1)+cn,dof*(k-1)+cm) = H(dof*(j-1)+cn,dof*(k-1)+cm) + a*b/(wr^2-w^2+i*2*w*wr*cr);
            
              end
          end
       end 
    end
end
