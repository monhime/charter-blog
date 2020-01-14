N = 400; %r方向の格子数
M = 300; %z方向の格子数
RMAX = 30; %r座標の最大値
ZMAX = 30; %z座標の最大値
Conv = 1e-6; %収束判定誤差
dr = RMAX / N; %r方向の刻み幅
dz = ZMAX / M; %z方向の刻み幅
phi = zeros(N,M); %0埋め 
flag = -10*ones(N,M); %非境界部分の-10埋め

for i = 1:N
    for j = 1:M
        %格子点(i,j)の座標(r,z)
        r = i*dr;
        z = j*dz;
        %ディリクレ条件(1kV): 高圧電極
        if ((r<1 && z > 16) || r^2+(z-14)^2<9)
            flag(i,j)=-3;
            phi(i,j)=1.0;
        %ディリクレ条件(0kV): 低圧電極、接地面、r遠方境界
        elseif ((r<1 && z < 5) || r^2+(z-7)^2<9 || j==1 || i==N)
            flag(i,j)=-3;
        %z方向ノイマン条件
        elseif (j == M)
            flag(i,j)=-1;
        %r方向ノイマン条件
        elseif (i == 1)
            flag(i,j)=-2;
        end
    end
end

MaxPhi = 1.0;
MaxErr = 1e20;

loop = 0;

r=repmat((dr:dr:N*dr)',1,M);%点(r,z)のrの値の行列

while MaxErr > Conv
    Prev_phi=phi;
        
    %ずらした行列を生成
    phi_l=circshift(phi,1,1);%r方向に1ずらした行列
    phi_h=circshift(phi,-1,1);%-r方向に1ずらした行列
    phi_j=circshift(phi,-1,2);%-zに1ずらした行列
    phi_k=circshift(phi,1,2);%zに1ずらした行列

    %各条件での計算結果を予め作成
    %非境界点
    phi_1=0.25*((phi_l+phi_h+phi_j+phi_k)+0.5*dr*(phi_h-phi_l)./r);    
    %z方向ノイマン条件
    phi_2=0.25 *((phi_l+phi_h+2*phi_k)+0.5*dr*(phi_h-phi_l)./r);     
    %r方向ノイマン条件
    phi_3=0.25*(2*phi_h+phi_j+phi_k);
    
    %非境界部分の点：flag(r,z)==-10
    %z方向ノイマン条件：flag(r,z)==-1
    %r方向ノイマン条件：flag(r,z)==-2
    %ディリクレ条件 (flag(r,z)==-3)は前回と同じであることに注意
    phi=(flag==-10).*phi_1+(flag==-1).*phi_2+(flag==-2).*phi_3+(flag==-3).*Prev_phi;
    
    %前回の行列との誤差。MaxPhiで規格化
    Err = (abs(phi-Prev_phi))/MaxPhi;
    MaxErr=max(max(Err));%行列の最大の要素
    loop = loop + 1;
    %1000ループごとにMaxMaxErr表示
    if rem(loop,1000)==0
       disp(loop);
       disp(MaxErr);
    end
end

rdata=linspace(dr,dr*N,N);
zdata=linspace(dz,dz*M,M);
createfigure(rdata,zdata,phi');

function createfigure(xdata1,ydata1,zdata1)
    figure1 = figure; % figure を作成
    axes1 = axes('Parent',figure1);% axes を作成 
    hold(axes1,'on');
    % contour を作成
    [c1,h1] = contour(xdata1,ydata1,zdata1,'TextStep',0.1,'LevelStep',0.01);
    clabel(c1,h1);
    ylabel({'z方向'}); % ylabel を作成
    xlabel({'r方向'}); % xlabel を作成
    box(axes1,'on');
    axis(axes1,'tight');
    % 残りの座標軸プロパティの設定
    set(axes1,'BoxStyle','full','FontName','Times New Roman','FontSize',24,'Layer','top','XGrid','on','XMinorGrid','on','XTick',[5 10 15 20 25 30],'YGrid','on','YMinorGrid','on','YTick',[5 10 15 20 25 30]);
    colorbar(axes1); % colorbar を作成 
end