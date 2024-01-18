function [results, timecc, spect] = LPfun_noplots(data, varargin)

ncomponents = [];
if length(varargin) == 1
    ncomponents = varargin{1};
end

N=size(data,1);
x=data(:,2);
t=data(:,1)-data(1,1);

deltat=t(2)-t(1);

%set up matrix from data (N-M)xM. We take M=0.75*N. It is backward prediction
M=floor(0.75*N);
X=zeros(N-M,M);

for i=1:(N-M)
   X(i,:)=x(i+1:i+M);
end
   
%compute the SVD of X to find the LS fit for linear prediction
[U, D, V] = svd(X);
d=diag(D);

d_length=length(d);
if isempty(ncomponents)
    semilogy(d,'o');  
    ncomponents=input('number of singular values to include?');
end

% "inverse" of rectangular diagonal matrix
E = zeros(M, N-M); % should be of size = zeros(size(D'));

for j=1:ncomponents
  E(j,j)=1/d(j);
end

xvector=x(1:N-M);
A = V*E*U'*xvector;

% polynomial roots of Z*A where Z = (z^2K, a1*z^(2K-1), a2*z^(2K-2),..., 1)
%roots sorted in descending order
s=sort(roots([1,-A']),'descend');  % sorting may be not necessary

% roots in the form z = exp(-b + i*w*dt)
b=log(abs(s(1:d_length)))/deltat;
[w,I]=sort(angle(s(1:d_length))/deltat);
b=b(I);

% remove double solutions and solutions with w==0
Nzeros=sum(w==0);
n_select=round((length(w)-Nzeros)/2+Nzeros);

w = abs(w(1:n_select));
b = b(1:n_select);

% keep solutions outsize of unit circle (when predicting backwards)
bnonneg=(b>=0);
b = b(bnonneg);
w = w(bnonneg);

% sort b's in decending order
[b,I]=sort(b,'descend');
w = w(I);
nw = length(w);

% This is X in terms of b's and w's
Xbar=zeros(N,nw*2+1);
t = (1:N)'*deltat;
for j=1:nw
    Xbar(:,2*j-1)=exp(-b(j).*t).*cos(w(j).*t);
    Xbar(:,2*j)=-exp(-b(j).*t).*sin(w(j).*t);
end
Xbar(:,nw*2+1)=ones(N,1);

% this finds the amplitudes and phases of the decaying cosines from the
% data
AA=lsqlin(Xbar,x);
for i=1:nw
    if AA(2*i-1)==0 && AA(2*i)==0
        C(i)=0;
        fi(i)=0;
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

% assemble the output and prepare the plots
yy(1:N)=0;
reconst=yy';
time=zeros(N,3+nw);
time(:,1:2)=[t,x];
timecc=zeros(N,nw);
Wmax=max(w);
Wmax=Wmax/2/pi;
bi=b/2/pi;
ww=w/2/pi;
freq = linspace(0, 1.5*Wmax, 500);
freq=freq';
spect=zeros(length(freq),1+nw);
spect(:,1)=freq;
spect_1=zeros(length(freq),nw);

totals=0;

for i=1:nw
    timecc(:,i)=C(i).*exp(-b(i).*t).*cos(w(i).*t+fi(i));
    reconst=reconst+timecc(:,i);
    time(:,3+i)=timecc(:,i);
    spect(:,1+i)=C(i)*bi(i)./((ww(i)-freq).^2+bi(i)^2);
    totals=totals+spect(:,1+i);
    spect_1(:,i)=spect(:,1+i);
end
time(:,3)=reconst;
Chi2=sum((reconst-x).^2)/N
reconst=reconst+AA(2*nw+1);

%% create a big figure with 6 panels
plotM=2; 
plotN=3;
subplot(plotN,plotM,1)
plot(t,timecc,t,x);
title('components', 'FontWeight','normal');
subplot(plotN,plotM,2)
plot(t,reconst,t,x);
title('LP fit + data','FontWeight','normal');
subplot(plotN,plotM,3);
% xl=s3.XLim;
% yl=s3.YLim;
% s3.XLimMode='manual';
% s3.YLimMode='manual';
plot(freq,spect_1,freq,totals);
% s3.XLim = xl;
% s3.YLim = yl;
title('LP spectrum', 'FontWeight','normal');

W=w'/2/pi;
B=b'/2/pi;
results=[W 
   B
   fi
   C];

results=results';
disp('        W(THz)      B          Phi       C')
disp(results)
    
residue=x-reconst;
subplot(plotN,plotM,4)
plot(t,residue);
title('Residue', 'FontWeight','normal');
xlabel('time (ps)');
ylabel('Data-Fit');

%FFt of the residue
[ws, Xw] = tdsfft(t, residue);
subplot(plotN,plotM,5);
% xl=s5.XLim;
% s5.XLimMode='manual';
% axis manual
plot(ws,abs(Xw));
% s5.XLim = xl;
title('FFT of the Residue', 'FontWeight','normal');
xlabel('freq (THz)');
ylabel('FFT');
