%1ά2���ۺϳ���
clc;
clear;
close all;

%Ƶ�ʷ�Χ
fre = 500;
ww1 = 2 * pi * 1;
ww2 =  2 * pi * fre;

%���ɶ�
dof = 1;
%%
%��������
%����ĸ���
manum=2;
%ÿ�������µ��ɵĸ���
nsu=zeros(manum,1);
for qq=1:manum
     nsu(qq,1)=4*dof;
end
%�²��ӽṹ֮��ĵ������ӵ���Ŀ
shnum1=4*dof;
%�ϲ��ӽṹ֮��ĵ������ӵ���Ŀ
shnum2 = 8*dof;

%�ӽṹA���ڵ���
inumA = 8*dof;
%�ӽṹB���ڵ���
inumB = 0*dof;
%�ӽṹc���ڵ���
inumC = 2*dof;

%ģ̬����
cr = 0;
%%

HCCABs = zeros(shnum1,shnum1,fre);
HCCAs = zeros(shnum1,shnum1,fre);
HICABs = zeros(inumA,shnum1,fre);
HIIAs =zeros(inumA,inumA,fre);
HICAAs = zeros(inumA,shnum1,fre);

HccBs = zeros(shnum1,shnum1,fre);
HddAs = zeros(shnum1,shnum1,fre);
HicAs = zeros(inumA,shnum1,fre);
HiiAs =zeros(inumA,inumA,fre);

kz1z = ones(1,dof);
kzy1z = ones(1,dof);
kz2z = ones(1,dof);
kzy2z = ones(1,dof);
%%
%�²����������
kz1z = 1610000*kz1z;
skz1z = [kz1z kz1z kz1z kz1z]; 
kzy1z = 0*kz1z;
skzy1z = [kzy1z kzy1z kzy1z kzy1z];
Z111 = diag(skz1z.*(1+i*skzy1z));%ע���ǳ���-w*w��ע����Щ�����ĵ�λ��g/N
Z112 = - diag(skz1z.*(1+i*skzy1z));
Z122 = diag(skz1z.*(1+i*skzy1z));

Z_UpUp1 = Z111;
Z_UtoD1 = Z112;
Z_DtoU1 = Z112;
Z_DownD1 = Z122; 
%%
%�ϲ�������迹
kz2z = 483000*kz2z;
skz2z = [kz2z kz2z kz2z kz2z kz2z kz2z kz2z kz2z]; 
kzy2z = 0*kzy2z;
skzy2z = [kzy2z kzy2z kzy2z kzy2z kzy2z kzy2z kzy2z kzy2z];

Z211 = diag(skz2z.*(1+i*skzy2z));%ע���ǳ���-w*w��ע����Щ�����ĵ�λ��g/N
Z212 = - diag(skz2z.*(1+i*skzy2z));
Z222 = diag(skz2z.*(1+i*skzy2z));

Z_UpUp2 = Z211;
Z_UtoD2 = Z212;
Z_DtoU2 = Z212;
Z_DownD2 = Z211; 
%%
%��������ӽṹ��ģ̬��Ϣ
% load('modaldata.mat');
%�ӽṹA��Ƶ�������ļ���������
fA=load('fA');%�������ļ�������ansysģ̬����֮������ļ�ֱ������
mA=load('subA.txt');%�������ӵ��ģ̬λ��
%�ӽṹB��Ƶ�������ļ�(����)
fB=load('fB');%�������ļ�������ansysģ̬����֮������ļ�ֱ������
mB=load('subB.txt');%������6�����ӵ��ģ̬λ�ƣ��ڲ������������ʱ��ֻҪ6�����ӵ��ģ̬λ�������㴫�ݺ���������Ҫ�����ڲ����ģ̬λ��
%�ӽṹC��Ƶ�������ļ������飩
fC=load('fC');%�������ļ�������ansysģ̬����֮������ļ�ֱ������
mC=load('subC.txt');%������6�����ӵ��ģ̬λ�ƣ��ڲ������������ʱ��ֻҪ6�����ӵ��ģ̬λ�������㴫�ݺ���������Ҫ�����ڲ����ģ̬λ��
%%
flag = 1;

for ww = ww1 :  2 * pi : ww2    
    wws(flag,1) = ww;
%�ӽṹ��Ƶ�����
%��ansys�õ��������ģ̬���ݣ���ͨ��ģ̬���ӵõ�Ƶ�캯��
%%
%�ӽṹB��4�����ӵ��Ƶ�죨������
HccB = calhFLEXIBLE(fB,mB,mB,shnum1,shnum1,0,0,ww,cr,dof);
HccBs(:,:,flag) = HccB;    
%%
%�ӽṹA���ڵ�֮���Ƶ�죨������    
HiiA = calhFLEXIBLE(fA,mA,mA,inumA,inumA,0,0,ww,cr,dof);
HiiAs(:,:,flag) = HiiA; 
%�ӽṹA�ڵ�����ӵ��Ƶ��   
HicA = calhFLEXIBLE(fA,mA,mA,inumA,shnum1,0,inumA,ww,cr,dof);
HicAs(:,:,flag) = HicA;
%�ӽṹA�����ӵ�֮���Ƶ��
HccA = calhFLEXIBLE(fA,mA,mA,shnum1,shnum1,inumA,inumA,ww,cr,dof);
HccAs(:,:,flag) = HccA;
%�ӽṹA���ӵ���ڵ��Ƶ��    
HciA = conj(HicA');
%%
%�ӽṹD�����飩
HiiD(1:inumC,1:inumC) = calhFLEXIBLE(fC,mC,mC,inumC,inumC,0,0,ww,cr,dof);
HiiD(inumC+1:2*inumC,inumC+1:2*inumC) = calhFLEXIBLE(fC,mC,mC,inumC,inumC,0,0,ww,cr,dof);

%������Ըĳ�ѭ��
HicD(1:inumC,1:nsu(1,1)) =calhFLEXIBLE(fC,mC,mC,inumC,nsu(1,1),0,inumC,ww,cr,dof);
HicD(inumC+1:2*inumC,nsu(1,1) + 1 : 2*nsu(1,1)) =calhFLEXIBLE(fC,mC,mC,inumC,nsu(1,1),0,inumC,ww,cr,dof);      
HciD = conj(HicD');

HccD(1:nsu(1,1),1:nsu(1,1)) = calhFLEXIBLE(fC,mC,mC,nsu(1,1),nsu(1,1),inumC,inumC,ww,cr,dof);
HccD(nsu(1,1) + 1 : 2*nsu(1,1) ,nsu(1,1) + 1 : 2*nsu(1,1)) = calhFLEXIBLE(fC,mC,mC,nsu(1,1),nsu(1,1),inumC,inumC,ww,cr,dof);
%%  
%��һ���ۺ�
PZUD1 = inv(Z_UtoD1);
PZDU1 = inv(Z_DtoU1);
     
D1 = HccB * (Z_DownD1 * PZUD1 * Z_UpUp1 - Z_DtoU1) * HccA + (HccB * Z_DownD1 * PZUD1 + PZUD1 * Z_UpUp1 * HccA ) + PZUD1;
PD1 = inv(D1);
      
HIIA = HiiA + HicA * PD1 * (HccB * Z_DtoU1 - HccB * Z_DownD1 * PZUD1 * Z_UpUp1 -  PZUD1 * Z_UpUp1) * HciA;
HICAA = HicA * PD1 * (HccB * Z_DownD1 + eye(size(HccB * Z_DownD1))) * PZUD1;
HICAB = - HicA * PD1 * HccB;
HCCA = HccA * PD1 * (HccB * Z_DownD1 + eye(size(HccB * Z_DownD1))) * PZUD1;
HCCAB = - HccA * PD1 * HccB;
HCCB = HccB * conj(PD1') * (HccA * Z_UpUp1 + eye(size(HccA * Z_UpUp1))) * PZDU1;

HCCABs(:,:,flag) = HCCAB;
HCCAs(:,:,flag) = HCCA;
HICABs(:,:,flag) = HICAB;
HIIAs(:,:,flag) = HIIA;
HICAAs(:,:,flag) = HICAA;

%��Ӧ����������ЩƵ�캯��������
HiiAB = [HCCA HCCAB;conj(HCCAB') HCCB];
HicAB = [conj(HICAA');conj(HICAB')];
HciAB = [HICAA HICAB];
HccAB = HIIA;
 
%%
%�ڶ����ۺϣ��ϲ��������
PZUD2 = inv(Z_UtoD2);
PZDU2 = inv(Z_DtoU2);
     
D2 = HccAB*(Z_DownD2 * PZUD2 * Z_UpUp2 - Z_DtoU2) * HccD + (HccAB * Z_DownD2 * PZUD2 + PZUD2 * Z_UpUp2 * HccD) + PZUD2;
PD2 = inv(D2);
     
HIIA2 = HiiD + HicD * PD2 * (HccAB * Z_DtoU2 - HccAB * Z_DownD2 * PZUD2 * Z_UpUp2 -  PZUD2 * Z_UpUp2) * HciD;
HIIAB2 = -HicD * PD2 * HciAB;
HICAA2 = HicD * PD2 * (HccAB * Z_DownD2 + eye(size(HccAB * Z_DownD2))) * PZUD2;
HICAB2 = - HicD * PD2 * HccAB;
     
HIIB2 = HiiAB + HicAB * conj(PD2') * (HccD * Z_UtoD2 - HccD * Z_UpUp2 * PZDU2 * Z_DownD2 -  PZDU2 * Z_DownD2) * HciAB;
HICBA2 = - HicAB * conj(PD2') * HccD;
HICBB2 = HicAB * conj(PD2') * (HccD * Z_UpUp2 + eye(size(HccD * Z_UpUp2))) * PZDU2;
     
HCCA2 = HccD * PD2 * (HccAB * Z_DownD2 + eye(size(HccAB * Z_DownD2))) * PZUD2;
HCCAB2 = - HccD * PD2 * HccAB;
     
HCCB2 = HccAB * conj(PD2') * (HccD * Z_UpUp2 + eye(size((HccD * Z_UpUp2)))) * PZDU2; 

%��Ӧ����������ЩƵ�캯��������
TotalH2 = [HIIA2 HIIAB2 HICAA2 HICAB2; conj(HIIAB2') HIIB2 HICBA2 HICBB2; conj(HICAA2') conj(HICBA2') HCCA2 HCCAB2; conj(HICAB2') conj(HICBB2') conj(HCCAB2') HCCB2];
HIIA2s(:,:,flag) = HIIA2;
HIIAB2s(:,:,flag) = HIIAB2;
HICAA2s(:,:,flag) = HICAA2;
HICAB2s(:,:,flag) = HICAB2;
%%
flag = flag + 1;
end
figure
ansHCCAB = load('HIIAB2.txt');
plot(1:500,log10(abs(ansHCCAB(:,2))));
hold on
syn(:,1) = HIIAB2s(1,8,:);
plot(log10(abs(syn(:,1))),':r','linewidth',2)
hold on