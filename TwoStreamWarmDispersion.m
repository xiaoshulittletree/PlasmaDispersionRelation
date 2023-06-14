##This scripts solved the equation (1) in [Shu Zhang et al., ion and electron bursts ... Nature Physics 2023]
##to show the calculation of growth rate's dispersion relation of IAI. 

Z=@(xi) (1i*sqrt(pi)*exp(-xi.^2).*(1+erf(sym(1i*xi))));
dZdx=@(xi) -2*(1+Z(xi).*xi);
%dZdx=@(xi) zprimeStevenc(xi);
wavenumber=0.1:0.01:0.8; %in the unit of 1/lambda_D, definition of LambdaD=sqrt(Te/me)/wpe
v0=0.12;%in the unit of sqrt(2Te/me) %0.12 for paper
Te=0.2;%keV %$0.2 for paper
Ti1=0.4; %keV %0.4 for paper
%Ti1=360;
Ti2=0.7;
vi1=0;
vi2=2;#for a secondary ion component
%R=1/100;
%R=1/1836;
R=1/1836/65; %for paper
%R=1/25;
Z=18;%for paper
%Z=1
n1=1; #density for first component of the ions
%Z=18;
wbunemanreal=zeros(size(wavenumber)); %in the unit of wpe?
wbunemanimag=wbunemanreal;
S(1)=0.0001i+0.001;
%S(1)=0.014;
for k=1:length(wavenumber)
    syms w;
    eqn=1-1/wavenumber(k)^2/2*(dZdx(w/wavenumber(k)/sqrt(2)-v0))...
       -n1/wavenumber(k)^2/2*Z*Te/Ti1*dZdx(w/wavenumber(k)/sqrt(2)*sqrt(Te/Ti1/R)-vi1)...
       -(1-n1)/wavenumber(k)^2/2*Z*Te/Ti2*dZdx(w/wavenumber(k)/sqrt(2)*sqrt(Te/Ti2/R)-vi2);
    %S=fsolve(eqn, 0.001+0.001i)
    if ~isempty(S)
        Sprevious=S(1);
        S=vpasolve(eqn,w,S(1))
    else
        %S=vpasolve(eqn,w,[1e-5*(1+1i), 1e-2*(1+1i)],'Random',false)
        S=vpasolve(eqn,w,Sprevious)

        disp(k)
        disp(wavenumber(k))
    end
    if ~isempty(S)
    wbunemanreal(k)=double(real(S(1)));
    wbunemanimag(k)=double(imag(S(1)));
    else
        wbunemanreal(k)=nan;
        wbunemanimag(k)=nan;
    end
end
figure;plot(wavenumber,wbunemanreal,wavenumber,wbunemanimag)

