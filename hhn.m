tend = 10;
dt = 0.01;
time = 0:dt:tend;
v1 = zeros(length(time),1);
voltage = -100:0.01:100;

n = zeros(length(time),1);
m = zeros(length(time),1);
h = zeros(length(time),1);
r = zeros(length(time),1);
n(1) = 0.0;
h(1) = 1.0;
m(1) = 0.0;
r(1) = 1.0
nn = zeros(length(voltage),1);
mm = zeros(length(voltage),1);
hh = zeros(length(voltage),1);
rr = zeros(length(voltage),1);

for k= 1:length(voltage)
    v = voltage(k);
    for i=2:length(time)
        n(i) = n(i-1) + (alphan(v)*(1-n(i-1)) - betan(v)*n(i-1))*dt;
        m(i) = m(i-1) + (alpham(v)*(1-m(i-1)) - betam(v)*m(i-1))*dt;
        h(i) = h(i-1) + (alphah(v)*(1-h(i-1)) - betah(v)*h(i-1))*dt;
        r(i) = r(i-1) + (rinf(v) - r(i-1))/tauh(v)*dt;
    end
nn(k) = n(length(time));
mm(k) = m(length(time));
hh(k) = h(length(time));
rr(k) = r(length(time));
end

%plot(voltage,nn,'LineWidth',3)
%hold on;
%plot(voltage,mm,'LineWidth',3)
%hold on;
%plot(voltage,hh,'LineWidth',3)
%hold on;
plot(voltage,rr,'LineWidth',3,'color',[0 0.5 0])
%legend('n','m','h','r')
xlabel('voltage (mV)')
ylabel('r')
set(gca,'fontsize',14)

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
function y = rinf(v)
    y = 1.0/(1.0+ exp((v+80.0)/10.2));
end

function taur = tauh(v)
    taur = 1.0/(exp(-14.59 - 0.086*v) + exp(-1.87+0.0701*v));
end
