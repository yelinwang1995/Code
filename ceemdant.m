%Function for  CEEMDANT

%Author: Dr. Yanpeng Cai, Dr. Youjie Li, Mr. Yelin Wang
%E-mail:yanpeng.cai@gdut.edu.cn(YP Cai), liyoujie@kust.edu.cn(YJ, Li), wangyelin0@163.com(YL Wang) 
%Last version: 30 aug 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CEEMDANLT is an optimized vision of iceemdan
%This code is developed based on the code of iceemdan.
%The file of iceemdan is available at http://www.bioingenieria.edu.ar/grupos/ldnlys/
%Colominas MA, Schlotthauer G, Torres ME. "Improve complete ensemble EMD: A suitable tool for biomedical signal processing" Biomedical Signal Processing and Control vol. 14 pp. 19-29 (2014)

%   Syntax
% modes=ceemdant(x,Nstd,NR,MaxIter)
% [modes its]=ceemdanlt(x,Nstd,NR,MaxIter)

%   Description

% OUTPUT
% modes: contain the obtained modes in a matrix with the rows being the modes        
% its: contain the sifting iterations needed for each mode for each realization (one row for each realization)

%INPUT
%x: signal to decompose
%Nstd: noise standard deviation
%NR: number of realizations
%MaxIter: maximum number of sifting iterations allowed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [modes,its]=ceemdanlt(x,Nstd,NR,MaxIter)
x=x(:)';
desvio_x=std(x);
x=x/desvio_x;
n=size(x,2);
modes=zeros(size(x));
temp=zeros(size(x));
aux=zeros(size(x));
iter=zeros(NR,round(log2(length(x))+5));

for i=1:NR
    white_noise{i}=randn(size(x));%creates the noise realizations
end;

for i=1:NR
    modes_white_noise{i}=emd(white_noise{i});%calculates the modes of white gaussian noise
end;

for i=1:NR %calculates the first mode
    xi=x+Nstd*modes_white_noise{i}(1,:)/std(modes_white_noise{i}(1,:));
    [temp, o, it]=emd(xi,'MAXMODES',1,'MAXITERATIONS',MaxIter);
    temp=temp(1,:);
    aux=aux+(xi-temp)/NR;
    iter(i,1)=it;
end;

%saves the temporary mode
modes_tem= x-aux; 
%%%limte the first mode by threshold
%calculate the threshold by quantile
q=quantile(modes_tem,3);
th_down=q(1)-1.5*(q(3)-q(1));
th_up=q(3)+1.5*(q(3)-q(1));
%de-nosing using the idea of threshold
for j=1:n
    if modes_tem(j)>=0
       tmp = (modes_tem(j)-th_up);
       modes_tem1(j)=(tmp+abs(tmp))/2;
    else
       tmp = (modes_tem(j)-th_down);
       modes_tem1(j)=(tmp-abs(tmp))/2;
    end
end
modes=modes_tem-modes_tem1;
medias = aux+modes_tem1;
 
k=1;
aux=zeros(size(x));
aux_remove=zeros(size(x));
es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));

while es_imf>1 %calculates the rest of the modes
    for i=1:NR
        tamanio=size(modes_white_noise{i});
        if tamanio(1)>=k+1
            noise=modes_white_noise{i}(k+1,:);
            noise=Nstd*noise;
            try
                [temp,o,it]=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch    
                it=0; disp('catch 1 '); disp(num2str(k))
                temp=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            end;
            temp=temp(end,:);
        else
            try
                [temp, o, it]=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch
                temp=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
                it=0; disp('catch 2 sin ruido')
            end;
            temp=temp(end,:);
        end;
        aux=aux+temp/NR;
    iter(i,k+1)=it;    
    end;
    modes=[modes;medias(end,:)-aux];
    medias = [medias;aux];
    aux=zeros(size(x));
    k=k+1;
    es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));
end;
modes = [modes;medias(end,:)];
modes=modes*desvio_x;
its=iter;
