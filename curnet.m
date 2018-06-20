%% 
%diffusion alone

nstep=100000;
dt=0.01;
T=nstep*dt;


numsynapse=6;
compstem=50;
numcomp=numsynapse*compstem;
V=zeros(numcomp,nstep);
I=1.;
gamma=10;
gld=0.001;


for j=1:nstep-1
    for i=1:numcomp
        if (i==1)
            V(i,j+1)=V(i,j)+(gamma*(V(i+1,j)-V(i,j))-gld*V(i,j)+I)*dt;
        else
            if (i==numcomp)
                V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)-V(i,j))-gld*V(i,j))*dt;
            else
                V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)+V(i+1,j)-2*V(i,j))-gld*V(i,j))*dt;
            end
        end
    end
end

plot(V')
mesh(V)

%%
%fluctuating inputs

rate=0.01;
taui=20*dt;
input=zeros(numsynapse,nstep);

%compute spike traces
for i=1:numsynapse
    itimes=zeros(1,nstep);
    time=0;
    while (time<T)
        time=time-log(rand())/rate;
        itimes(floor(time/dt))=1/dt;
    end
    for j=1:nstep
        input(i,j+1)=input(i,j) + (-input(i,j)/taui + itimes(j))*dt;
    end
end

%figure(2)
%plot(input')

%%
%diffusion with fluctuating inputs

for j=1:nstep-1
    for i=1:numcomp
        if (i==1)
            V(i,j+1)=V(i,j)+(gamma*(V(i+1,j)-V(i,j))-gld*V(i,j))*dt;
        else
            if (i==numcomp)
                V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)-V(i,j))-gld*V(i,j))*dt;
            else
                if (mod(i,compstem)==1)
                     V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)+V(i+1,j)-2*V(i,j))-gld*V(i,j)+input(floor(i/compstem),j))*dt;
                else
                     V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)+V(i+1,j)-2*V(i,j))-gld*V(i,j))*dt;
                end
            end
        end
    end
end

mesh(V)

%%
%diffusion with fluctuating conductance

gls=1;
Ers=2*double(randn(1,numsynapse)>0)-1;

for j=1:nstep-1
    for i=1:numcomp
        if (i==1)
            V(i,j+1)=V(i,j)+(gamma*(V(i+1,j)-V(i,j))-gld*V(i,j))*dt;
        else
            if (i==numcomp)
                V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)-V(i,j))-gld*V(i,j))*dt;
            else
                if (mod(i,compstem)==1)
                     isyn=floor(i/compstem);
                     V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)+V(i+1,j)-2*V(i,j))-gld*V(i,j)+gls*input(floor(i/compstem),j)*(Ers(isyn)-V(i,j)))*dt;
                else
                     V(i,j+1)=V(i,j)+(gamma*(V(i-1,j)+V(i+1,j)-2*V(i,j))-gld*V(i,j))*dt;
                end
            end
        end
    end
end

mesh(V)
