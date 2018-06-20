tend = 1000;
dt = 0.01;
time = 0:dt:tend;
trials = 10;
gsyn1 = 4.0;
gsyn2 = 2.0;
esyn2 = -70.0;
esyn1 = 0.0;
gama = 1.0;
ena = 50.0;
ek = -77.0;
el = -54.4;
gna = 120.0;
gk = 36.0;
gl = 0.3;
sigma2 =0.1:0.1:5.0;
vmean = zeros(length(sigma2),trials);
I =10.0;
C = 0.001;
taus = 2.0;
for k= 1:length(sigma2)
    sigma = sigma2(k)
    for g=1:trials
    spikes = 0;
    t = 0;
    v2 = zeros(length(time),1);
n = zeros(length(time),1);
m = zeros(length(time),1);
h = zeros(length(time),1);
n(1) = 0.0;
h(1) = 1.0;
m(1) = 0.0;
v2(1) = .9;

w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
frate = 1.0 ;
    sp1 = zeros(length(time),1);
    sp2 = zeros(length(time),1);
    s1 = zeros(length(time),1);
    s2 = zeros(length(time),1);
    times = zeros(length(time),1);
    while t < tend
    %t2=t
        tw = -log(rand())/frate;
        t = t + tw;
        t1 = t;
        if t1 < 0
            t1 = 0.0;
        end
        
        t2 = t + sigma*randn();
        if t2< 0
            t2 = 0.0;
        end
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1) + sp1(i-1) )*dt/taus ;
        s2(i) = s2(i-1) + (-s2(i-1) + sp2(i-1) )*dt/taus ;
    end



    for i=2:length(time)
        n(i) = n(i-1) + (alphan(v2(i-1))*(1-n(i-1)) - betan(v2(i-1))*n(i-1))*dt;
        m(i) = m(i-1) + (alpham(v2(i-1))*(1-m(i-1)) - betam(v2(i-1))*m(i-1))*dt;
        h(i) = h(i-1) + (alphah(v2(i-1))*(1-h(i-1)) - betah(v2(i-1))*h(i-1))*dt;
        v2(i) = v2(i-1)+  ( (gna*(m(i-1)^3)*h(i-1)*(-v2(i-1)+ena)) + (gk*(n(i-1)^4)*(-v2(i-1)+ ek)) + (gl*(-v2(i-1) + eleak))  + ( gama*(w2(i-1) - v2(i-1))))*dt;
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i-1)) + (gl*(-v1(i-1) + eleak)))*dt;
        w1(i) = w1(i-1) + (gsyn1*s1(i-1)*(esyn1 - w1(i-1)) +  gama*(v1(i-1) - w1(i-1)) + (gl*(-w1(i-1) + eleak)))*dt;
        w2(i) = w2(i-1) + (gsyn2*s2(i-1)*(esyn2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ (gl*(-w2(i-1) + eleak)) )*dt;
    end
    
for l=2:(length(time)-1)
        if v2(l-1)< v2(l)&& v2(l) > v2(l+1)&&v2(l)>0.50
            spikes = spikes + 1;
            times(l) = 1;
        end
    end
vmean(k,g) = spikes;
    end
end
figure(1)
plot(time,v2);
hold on;
plot(time,v1);
%plot(time,s2)
ylabel('voltage (mV)')
xlabel('time (msec)')
set(gca,'fontsize',14)

figure(2)
plot(sigma2,mean(vmean,2),'LineWidth',3)
xlabel('jitter')
ylabel('Spiking rate')
function y = alpham(v) 
    y = 0.1*(v+40.0)/(1.0 - exp(-0.1*(v+40.0)));
end
function y = betam(v)
    y = 4.0*exp(-(v+65.0)/18.0);
end
function y = alphan(v)
    y = 0.01*(v+55.0)/(1.0 - exp(-(v+55.0)/10.0));
end
function y = betan(v)
    y = 0.125*exp(-(v+65.0)/80.0);
end
function y = alphah(v)
    y = 0.07*exp(-(v+65.0)/20.0);
end
function y = betah(v)
    y = 1.0/(1.0+ exp(-(v+35.0)/10.0));
end
