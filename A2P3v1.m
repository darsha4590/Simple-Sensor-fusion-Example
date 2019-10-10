clc;
clear all;

X=linspace(-100,100,401);
Y=linspace(-100,100,401);
XY=[];
for i=1:length(X)
    %dum=0;
    for j=1:length(Y)
        XY=[XY;X(1,i) Y(1,j)];
    end
end

%% Constants and given parameter values
th1=0.453;
th2=2.3238;
sig1=0.1;
sig2=0.2;
S1=[1 1];
S2=[4,1];

%% We will be evaluating the values for MAP estimate
Z=[];
for i =1:length(XY(:,1))
    x=XY(i,1);
    y=XY(i,2);
    Ph1=(th1-atan2((y-S1(1,2)),(x-S1(1,1))));
    Ph2=(th2-atan2((y-S2(1,2)),(x-S2(1,1))));
    Map=exp(-1*(((Ph1^2/sig1^2)+(Ph2^2/sig2^2))));
    %Map=(((Ph1^2/sig1^2)+(Ph2^2/sig2^2)));
    Z=[Z; Map];
end
plot(1)
scatter3(XY(:,1),XY(:,2),Z)
[k,id]=max(Z);
%Map estimate value in variable below
Mapestimate=XY(id,1:2)

%[Xm Ym]=ind2sub(size(Z),k)
tanf1=@(x,y)atan2((y-S1(1,2)),(x-S1(1,1)))
tanf2=@(x,y)atan2((y-S2(1,2)),(x-S2(1,1)))
%first distribution function defined
fun1=@(x,y)x*((2*pi*sig1*sig2)^-1)*exp(-1*(((th1-tanf1(x,y))^2/sig1^2)+((th2-tanf2(x,y))^2/sig2^2)));
%second distribution function defined
fun2=@(x,y)((2*pi*sig1*sig2)^-1)*exp(-1*(((th1-tanf1(x,y))^2/sig1^2)+((th2-tanf2(x,y))^2/sig2^2)));
fun3=@(x,y)y*((2*pi*sig1*sig2)^-1)*exp(-1*(((th1-tanf1(x,y))^2/sig1^2)+((th2-tanf2(x,y))^2/sig2^2)));
%fun1=@(x,y)1
%q=integral2(fun1,-100,100,-100,100)
ezsurf(fun2);

%%numerical integration for the bayes estimate
by1=linspace(-100,100,401);
bx1=linspace(-100,100,401);
dx=bx1(1,2)-bx1(1,1);
dy=by1(1,2)-by1(1,1);
dum=[];
%numerical integration is done using flat integration scheme
for i=1:length(by1)-1
    for j=1:length(bx1)-1
        %Bayes estimate for x
        z=0.25*(fun1(bx1(1,i),by1(1,j))+fun1(bx1(1,i+1),by1(1,j))+fun1(bx1(1,i),by1(1,j+1))+fun1(bx1(1,i+1),by1(1,j+1)));
        %this is used for the normalizing factor evaluation
        z2=0.25*(fun2(bx1(1,i),by1(1,j))+fun2(bx1(1,i+1),by1(1,j))+fun2(bx1(1,i),by1(1,j+1))+fun2(bx1(1,i+1),by1(1,j+1)));
        %Bayes estimate for y
        z3=0.25*(fun3(bx1(1,i),by1(1,j))+fun3(bx1(1,i+1),by1(1,j))+fun3(bx1(1,i),by1(1,j+1))+fun3(bx1(1,i+1),by1(1,j+1)));
        dum=[dum;(z*dx*dy) (z2*dx*dy) (z3*dx*dy)];
    end
end
%Bayes estimation value below
Bayesestimate=[(sum(dum(:,1))/sum(dum(:,2))) (sum(dum(:,3))/sum(dum(:,2)))];
%m=[bx1(1,124001) by1(1,124001)];

        
















