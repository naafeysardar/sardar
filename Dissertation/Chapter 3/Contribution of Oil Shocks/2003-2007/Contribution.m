clear all;

global t

load oil.txt; o=oil; [t,q]=size(o);

load oilsup.txt; x=oilsup; [t,q]=size(x);
load oilad.txt; y=oilad; [t,q]=size(y);
load oildd.txt; z=oildd; [t,q]=size(z);

p=0

irfsup = x(:,1)
oilsup = x(:,2)
irfsup = transpose(irfsup)
oilsup = transpose(oilsup)

irfad = y(:,1)
oilad = y(:,2)
irfad = transpose(irfad)
oilad = transpose(oilad)

irfdd = z(:,1)
oildd = z(:,2)
irfdd = transpose(irfdd)
oildd = transpose(oildd)

time=(2003+1/4:1/4:2007+4/4)';

xhat=zeros(t-p,1); yhat=zeros(t-p,1); zhat=zeros(t-p,1);
for i=1:t-p
    xhat(i,:)=dot(irfsup(1,1:i),oilsup(1,i:-1:1));
    yhat(i,:)=dot(irfad(1,1:i),oilad(1,i:-1:1));
    zhat(i,:)=dot(irfdd(1,1:i),oildd(1,i:-1:1));
end;

% Contribution to real price of oil
time=2003+1/4:1/4:2007+4/4;

subplot(2,1,1)
plot(time,yhat,'b-');
title('Cumulative Effect of Emerging Market Demand Shock on Real Price of Crude Oil')
axis([2003+1/4 2007+4/4 0 +0.5])
grid on

subplot(2,1,2)
plot(time,oil,'b-');
title('Real Price of Crude Oil')
axis([2003+1/4 2007+4/4 10 +45])
grid on


