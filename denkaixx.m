N = 400; %r�����̊i�q��
M = 300; %z�����̊i�q��
RMAX = 30; %r���W�̍ő�l
ZMAX = 30; %z���W�̍ő�l
Conv = 1e-6; %��������덷
dr = RMAX / N; %r�����̍��ݕ�
dz = ZMAX / M; %z�����̍��ݕ�
phi = zeros(N,M); %0���� 
flag = -10*ones(N,M); %�񋫊E������-10����

for i = 1:N
    for j = 1:M
        %�i�q�_(i,j)�̍��W(r,z)
        r = i*dr;
        z = j*dz;
        %�f�B���N������(1kV): �����d��
        if ((r<1 && z > 16) || r^2+(z-14)^2<9)
            flag(i,j)=-3;
            phi(i,j)=1.0;
        %�f�B���N������(0kV): �ሳ�d�ɁA�ڒn�ʁAr�������E
        elseif ((r<1 && z < 5) || r^2+(z-7)^2<9 || j==1 || i==N)
            flag(i,j)=-3;
        %z�����m�C�}������
        elseif (j == M)
            flag(i,j)=-1;
        %r�����m�C�}������
        elseif (i == 1)
            flag(i,j)=-2;
        end
    end
end

MaxPhi = 1.0;
MaxErr = 1e20;

loop = 0;

r=repmat((dr:dr:N*dr)',1,M);%�_(r,z)��r�̒l�̍s��

while MaxErr > Conv
    Prev_phi=phi;
        
    %���炵���s��𐶐�
    phi_l=circshift(phi,1,1);%r������1���炵���s��
    phi_h=circshift(phi,-1,1);%-r������1���炵���s��
    phi_j=circshift(phi,-1,2);%-z��1���炵���s��
    phi_k=circshift(phi,1,2);%z��1���炵���s��

    %�e�����ł̌v�Z���ʂ�\�ߍ쐬
    %�񋫊E�_
    phi_1=0.25*((phi_l+phi_h+phi_j+phi_k)+0.5*dr*(phi_h-phi_l)./r);    
    %z�����m�C�}������
    phi_2=0.25 *((phi_l+phi_h+2*phi_k)+0.5*dr*(phi_h-phi_l)./r);     
    %r�����m�C�}������
    phi_3=0.25*(2*phi_h+phi_j+phi_k);
    
    %�񋫊E�����̓_�Fflag(r,z)==-10
    %z�����m�C�}�������Fflag(r,z)==-1
    %r�����m�C�}�������Fflag(r,z)==-2
    %�f�B���N������ (flag(r,z)==-3)�͑O��Ɠ����ł��邱�Ƃɒ���
    phi=(flag==-10).*phi_1+(flag==-1).*phi_2+(flag==-2).*phi_3+(flag==-3).*Prev_phi;
    
    %�O��̍s��Ƃ̌덷�BMaxPhi�ŋK�i��
    Err = (abs(phi-Prev_phi))/MaxPhi;
    MaxErr=max(max(Err));%�s��̍ő�̗v�f
    loop = loop + 1;
    %1000���[�v���Ƃ�MaxMaxErr�\��
    if rem(loop,1000)==0
       disp(loop);
       disp(MaxErr);
    end
end

rdata=linspace(dr,dr*N,N);
zdata=linspace(dz,dz*M,M);
createfigure(rdata,zdata,phi');

function createfigure(xdata1,ydata1,zdata1)
    figure1 = figure; % figure ���쐬
    axes1 = axes('Parent',figure1);% axes ���쐬 
    hold(axes1,'on');
    % contour ���쐬
    [c1,h1] = contour(xdata1,ydata1,zdata1,'TextStep',0.1,'LevelStep',0.01);
    clabel(c1,h1);
    ylabel({'z����'}); % ylabel ���쐬
    xlabel({'r����'}); % xlabel ���쐬
    box(axes1,'on');
    axis(axes1,'tight');
    % �c��̍��W���v���p�e�B�̐ݒ�
    set(axes1,'BoxStyle','full','FontName','Times New Roman','FontSize',24,'Layer','top','XGrid','on','XMinorGrid','on','XTick',[5 10 15 20 25 30],'YGrid','on','YMinorGrid','on','YTick',[5 10 15 20 25 30]);
    colorbar(axes1); % colorbar ���쐬 
end