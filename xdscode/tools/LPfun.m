function [results, timecc, zz] = LPfun(data, varargin)

ncomponents = [];
if length(varargin) == 1
    ncomponents = varargin{1};
end

N=size(data,1);

if  mod(N,2)~=0
    N=N-1;
end

x=data(1:N,2);
t=data(1:N,1)-data(1,1);

deltat=t(2)-t(1);

M=N/2;

%set up matrix from data (N-M)xM. We take M=0.5*N. It is backward prediction

X=zeros(M,M);

for i=1:M
   X(:,i)=x(i+1:i+M);
end
 
%computation of the (N-M)x(N-M) noonegativmatrix XX'and diagonalization

XX=X*X';
[U,D] = eig(XX);
d=diag(D);

d_length=length(d);
no=1:d_length; 
if isempty(ncomponents)
    semilogy(no,d,'o');  
    ncomponents=input('number of singular values to include?');
end

E=zeros(M); 

for j=d_length-ncomponents+1:d_length %for j=1:size(D,1)
  E(j,j)=1/d(j);% 1 overlambda
end



xvector=x(1:M);
A=X'*U*E*U'*xvector;



%polynomial roots

r=roots([1,-A']);
s=sort(r);  % sorting may be not necessary

BB=length(s);
         
  %roots sorted in descending order

        ss(1:BB)=s(BB:-1:1);

        ss=ss';

        b(1:d_length)=log(abs(ss(1:d_length)))/deltat;
        w(1:d_length)=angle(ss(1:d_length))/deltat;

        [P,I]=sort(w);
        Z=b(I);
       
        %count number of w==0
        Nzeros=sum(w==0);

        n_select=round((d_length-Nzeros)/2+Nzeros);

        WW(1:n_select)=abs(P(1:n_select));
        BBB(1:n_select)=(Z(1:n_select));  

     %counting for positive damping constants  

        Npos=sum(BBB>=0);
        [B1,J]=sort(BBB);
        W1=WW(J);
         
  % sorted in descending order

        B1_length=length(B1);
        B2(1:B1_length)=B1(B1_length:-1:1);

        W2(1:B1_length)=W1(B1_length:-1:1);

        W(1:Npos)=W2(1:Npos);
        B(1:Npos)=B2(1:Npos);

        W_length=length(W);

        Xbar=zeros(N,W_length*2+1);
        i=1:N;
        i=i';
        for j=1:W_length;
           Xbar(:,2*j-1)=exp(-B(j).*(i-1)*deltat).*cos(W(j).*(i-1)*deltat);
           Xbar(:,2*j)=-exp(-B(j).*(i-1)*deltat).*sin(W(j).*(i-1)*deltat);
             
        end
        Xbar(:,W_length*2+1)=ones(N,1);

         
         
           

AA=lsqlin(Xbar,x);

for i=1:W_length

     if AA(2*i-1)==0 & AA(2*i)==0
           C(i)=0;fi(i)=0;
     elseif AA(2*i-1)==0
            fi(i)=sign(AA(2*i))*pi/2;
             C(i)=abs(AA(2*i));
     elseif AA(2*i)==0
             fi(i)=(1-sign(AA(2*i-1)))*pi/2;
             C(i)=abs(AA(2*i-1));
        else
   fi(i)=atan2(AA(2*i),AA(2*i-1));
 
    C(i)=sqrt(AA(2*i)^2+AA(2*i-1)^2);
    end
end





yy(1:N)=0;
zz=yy';
time=zeros(N,3+W_length);
time(:,1:2)=[t,x];
%timecc=zeros(N,1);
timecc=zeros(N,W_length);
Wmax=max(W);
%Wmax=Wmax*100/2/3/pi; 
%bi=B*100/2/3/pi;
%ww=W*100/2/3/pi;
Wmax=Wmax/2/pi;
bi=B/2/pi;
ww=W/2/pi;
% freq=0:Wmax/1000:1.5*Wmax;%in cm-1
freq = linspace(0, 1.5*Wmax, 500);
freq=freq';
spect=zeros(length(freq),1+W_length);
spect(:,1)=freq;
spect_1=zeros(length(freq),W_length);

totals=0;

for i=1:W_length
   
   %timecc=C(i).*exp(-B(i).*t).*cos(W(i).*t+fi(i));
   timecc(:,i)=C(i).*exp(-B(i).*t).*cos(W(i).*t+fi(i));
   %zz=zz+C(i).*exp(-B(i).*t).*cos(W(i).*t+fi(i)); 
   zz=zz+timecc(:,i);
   time(:,3+i)=timecc(:,i);
   spect(:,1+i)=C(i)*bi(i)./((ww(i)-freq).^2+bi(i)^2);
   totals=totals+spect(:,1+i);
   spect_1(:,i)=spect(:,1+i);
end
time(:,3)=zz;
Chi2=sum((zz-x).^2)/N
zz=zz+AA(2*W_length+1);

%% create a big figure with 6 panels
plotM=2; 
plotN=3;
subplot(plotN,plotM,1)
plot(t,timecc,t,x);
title('components');
subplot(plotN,plotM,2)
plot(t,zz,t,x);
title('LP fit + data');
s3=subplot(plotN,plotM,3)
xl=s3.XLim;
yl=s3.YLim;
s3.XLimMode='manual';
s3.YLimMode='manual';
plot(s3,freq,spect_1,freq,totals);
s3.XLim = xl;
s3.YLim = yl;
title('LP spectrum');

W=W/2/pi;
B=B/2/pi;
results=[W 
   B
   fi
   C];

results=results';
disp(['        W(THz)      B          Phi       C'])
disp(results)
    
residue=x-zz;
subplot(plotN,plotM,4)
plot(t,residue);title('Residue');xlabel('time(ps)');ylabel('Data-Fit');

%FFt of the residue
[ws, Xw] = tdsfft(t, residue);
s5=subplot(plotN,plotM,5)
xl=s5.XLim;
s5.XLimMode='manual';
axis manual
plot(ws,abs(Xw));
s5.XLim = xl;
title('FFT of the Residue');
xlabel('freq (THz)');
ylabel('FFT');
