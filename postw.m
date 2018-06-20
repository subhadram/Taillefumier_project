ffrate = [0.01]%:0.1:1.0;
%t = 0;
dt = 0.1;
tend = 1000.0;
sigma = 0.01;
time = 0:dt:tend;
vthresh = 0.50;
vreset = -1
vspike = 2.0;
j = 2;
vmean = zeros(length(ffrate),1);
g_d = 1;
g_s1 = 0.50;
g_s2 = 0.025;
gama = 2.0;
g_l = 0.8;
e_d = 0.0;
e_r1 = 0.5;
e_r2 = -0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 30.0;

%t = 0
for k= 1:length(ffrate)
    spikes = 0;
    t = 0;
    frate = ffrate(k) ;
    sp1 = zeros(length(time),1);
sp2 = zeros(length(time),1);
s1 = zeros(length(time),1);
s2 = zeros(length(time),1);

while t < tend
    %t2=t
    tw = -log(rand())/frate;
    
    t = t + tw;
    t1 = round(t,0);
    t2 = t + sigma*randn();
    t2 = round(t2,0);
    j=2;
    sp1(floor(t1/dt)+1)=1/dt;
    sp2(floor(t2/dt)+1)=1/dt;
end
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end

t = 0;
end

v1 = zeros(length(time),1);
w1 = zeros(length(time),1);
v2 = zeros(length(time),1);
w2 = zeros(length(time),1);

for i = 2:length(time)
    v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_l*(-v1(i-1) + e_d))*dt ;
    v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_l*(-v2(i-1) + e_d))*dt;
    w1(i) = w1(i-1) + 0.001 +(0*g_s1*s1(i-1)*(e_r1 -w1(i-1)) + gama*(v1(i-1) - w1(i-1)) + g_l*(-w1(i-1) + e_d))*dt;
    %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
    w2(i) = w2(i-1) + 0.001 + (0*g_s2*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
end
%plot(time,s1)
plot(time,v1)