tend = 1000;
dt = 0.01;
time = 0:dt:tend;
v1 = zeros(length(time),1);
voltage = -80:0.1:80;
n = zeros(length(time),1);
m = zeros(length(time),1);
h = zeros(length(time),1);

nn = zeros(length(voltage),1);
mm = zeros(length(voltage),1);
hh = zeros(length(voltage),1);


for k= 1:length(voltage)
    v = voltage(k)
    for i=2:length(time)
        n(i) = n(i-1) + (alphan(v)*(1-n(i-1)) - betan(v)*n(i-1))*dt
        m(i) = h(i-1) + (alpham(v)*(1-m(i-1)) - betam(v)*m(i-1))*dt
        h(i) = h(i-1) + (alphah(v)*(1-h(i-1)) - betah(v)*h(i-1))*dt
    end
nn(k) = n(i)
mm(k) = m(i)
hh(k) = h(i)
end

plot(voltage,nn)
hold on;
plot(voltage,mm)
hold on;
plot(voltage,hh)
legend('n','m','h')

function y = alpham(v) 
    y = 0.1*(v+40)/(1 - exp(-0.1*(v+40)));
end

function y = alphan(v)
    y = 0.01*(v+55)/(1 - exp(-0.1*(v+55)));
end

function y = alphah(v)
    y = 0.07*exp((v+65)/20.0);
end

function y = betam(v)
    y = 4.0*exp(-(v+65)/18.0);
end

function y = betan(v)
    y = 0.125*exp(-(v+65)/80.0);
end

function y = betah(v)
    y = (v+35)/(1+ exp((v+35)/10));
end
