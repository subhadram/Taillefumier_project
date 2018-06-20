%%
%tonic inhibtion
ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;
tonic_inhi = -.01;
trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    for g=1:trials
    spikes = 0;
    t = 0;
    frate = ffrate(k); 
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
        t2 = t + sigma*randn();
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        
        w1(i) = w1(i-1) + (tonic_inhi +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
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
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')
hold on;

figure(2)
plot(time,v2)
hold on;
plot(time,s2)
xlabel('')
ylabel('Spiking rate')


%%
%correlated inhibition
ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;
tonic_inhi = 0.1;
trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
total_inhi = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    for g=1:trials
    spikes = 0;
    t = 0;
    frate = ffrate(k) ;
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
        t2 = t + sigma*randn();
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        x(i) = g_s1*s1(i-1)*(e_r1 - w1(i-1));
        w1(i) = w1(i-1) + (g_s1*s1(i-1)*(e_r1 - w1(i-1)) +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
    end
    
    for l=2:(length(time)-1)
        if v2(l-1)< v2(l)&& v2(l) > v2(l+1)&&v2(l)>0.50
            spikes = spikes + 1;
            times(l) = 1;
        end
    end
total_inhi(k,g) = mean(x);
vmean(k,g) = spikes;
    end
end
figure(1)
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')
hold on;


inhi = mean(total_inhi,2);


figure(3)
plot(time,v2)
hold on;
plot(time,s2)
%%
%control
ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;

trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    for g=1:trials;
    spikes = 0;
    t = 0;
    frate = ffrate(k) ;
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
        t2 = t + sigma*randn();
        
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        
        w1(i) = w1(i-1) + (0*g_s1*s1(i-1)*(e_r1 - w1(i-1)) +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
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
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')

figure(4)
plot(time,v2)
hold on;
plot(time,s2)
%%
%normalized tonic inhibition


ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;
trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    tonic_inhi = inhi(k)
    for g=1:trials
    spikes = 0;
    t = 0;
    frate = ffrate(k); 
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
        t2 = t + sigma*randn();
        
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        
        w1(i) = w1(i-1) + (tonic_inhi +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
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
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')
legend('Tonic inhibtion','Correlated inhibitory synapse', 'control','normalized tonic inhibtion')
hold on;

figure(5)
plot(time,v2)
hold on;
plot(time,s2)
xlabel('')
ylabel('Spiking rate')
%%
%uncorrelated inhibition
ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;
tonic_inhi = 0.1;
trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
%total_inhi = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    for g=1:trials
    spikes = 0;
    t = 0;
    tm = 0;
    frate = ffrate(k) ;
    
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
        t2 = t + sigma*randn();
        %t2 = round(t2,0);
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        %sp2(floor(t2/dt)+1)=1/dt;
    end
    while tm < tend
    %t2=t
        tw2 = -log(rand())/frate;
        tm = tm + tw2;
        t2 = tm;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    tm = 0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        %x(i) = g_s1*s1(i-1)*(e_r1 - w1(i-1));
        w1(i) = w1(i-1) + (g_s1*s1(i-1)*(e_r1 - w1(i-1)) +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
    end
    
    for l=2:(length(time)-1)
        if v2(l-1)< v2(l)&& v2(l) > v2(l+1)&&v2(l)>0.50
            spikes = spikes + 1;
            times(l) = 1;
        end
    end
%total_inhi(k,g) = mean(x);
vmean(k,g) = spikes;
    end
end
figure(1)
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')
hold on;


%inhi = mean(total_inhi,2);


figure(6)
plot(time,v2)
hold on;
plot(time,s2)
%% different rate 

ffrate = 0.01:0.01:0.2;
dt = 0.1;
tend = 1000.0;
sigma = 0.001;
time = 0:dt:tend;
vthresh = 0.0;
vreset = -1;
vspike = 1.0;
tonic_inhi = 0.1;
trials = 100;
j = 2;
vmean = zeros(length(ffrate),trials);
%total_inhi = zeros(length(ffrate),trials);
g_d = 0.4;
g_s1 = 0.40;
g_s2 = 0.45;
gama = 0.2;
e_d = 0.0;
e_r1 = -0.5;
e_r2 = 0.5;
e_l = 0.0;
taus = 20.0*dt;
tau = 1.0*dt;

%v1(1) = 0.001 ;

for k= 1:length(ffrate)
    for g=1:trials
    spikes = 0;
    t = 0;
    tm = 0;
    frate = ffrate(k) ;
    frate2 = ffrate(5);
    sp1 = zeros(length(time),1);
    sp2 = zeros(length(time),1);
    s1 = zeros(length(time),1);
    s2 = zeros(length(time),1);
    times = zeros(length(time),1);
    while t < tend
    %t2=t
        tw = -log(rand())/frate2;
        t = t + tw;
        t1 = round(t,1);
        
        j=2;
        sp1(floor(t1/dt)+1)=1/dt;
        %sp2(floor(t2/dt)+1)=1/dt;
    end
    while tm < tend
    %t2=t
        tw2 = -log(rand())/frate;
        tm = tm + tw2;
        t2 = tm;
        sp2(floor(t2/dt)+1)=1/dt;
    end
    t =0;
    tm = 0;
    for i=2:length(time)
        s1(i) = s1(i-1) + (-s1(i-1)/taus + sp1(i-1) )*dt ;
        s2(i) = s2(i-1) + (-s2(i-1)/taus + sp2(i-1) )*dt ;
    end
    v1 = zeros(length(time),1);
    w1 = zeros(length(time),1);
    v2 = zeros(length(time),1);
    w2 = zeros(length(time),1);
    
    

    for i = 2:length(time)
        v1(i) = v1(i-1) + (gama*(w2(i-1) + w1(i-1) - 2*v1(i)) + g_d*(-v1(i-1)))*dt;
        %v2(i) = v2(i-1) + (gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
        if v2(i-1)>=vspike
             v2(i)=v2(i-1) + ((vreset-v2(i-1))/tau)*dt;
         elseif v1(i-1)<vthresh
             v2(i) = v2(i-1) + ( gama*(w2(i-1) - v2(i-1)) + g_d*(-v2(i-1)))*dt;
         elseif v2(i-1)>= vthresh && v2(i-1)<vspike
           v2(i)=vspike;
            
%             spikes= spikes + 1;
          end
        %x(i) = g_s1*s1(i-1)*(e_r1 - w1(i-1));
        w1(i) = w1(i-1) + (g_s1*s1(i-1)*(e_r1 - w1(i-1)) +  gama*(v1(i-1) - w1(i-1)) + g_d*(-w1(i-1)))*dt;
        %w2(i) = w2(i-1) + (g_s2*s2(i-1)*(e_r2 - w2(i-1)) + gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_l*(-w2(i-1) + e_d) )*dt;
        w2(i) = w2(i-1) +  (g_s2*s2(i-1)*(e_r2 - w2(i-1)) +gama*(v1(i-1) + v2(i-1) - 2*w2(i-1))+ g_d*(-w2(i-1)))*dt;
    end
    
    for l=2:(length(time)-1)
        if v2(l-1)< v2(l)&& v2(l) > v2(l+1)&&v2(l)>0.50
            spikes = spikes + 1;
            times(l) = 1;
        end
    end
%total_inhi(k,g) = mean(x);
vmean(k,g) = spikes;
    end
end
figure(1)
plot(ffrate,mean(vmean,2),'LineWidth',3)
xlabel('Firing rate')
ylabel('Spiking rate')
legend('Tonic inhibtion','Correlated inhibitory synapse', 'control','normalized tonic inhibtion','uncorrelated same rate','uncorrelated fixed inhin rate')
hold on;


%inhi = mean(total_inhi,2);


figure(7)
plot(time,v2)
hold on;
plot(time,s2)