function [x1,y1,xe]=VIIM_initialdata(angle,N,eps,T)

angle=angle*pi/180;
coef=(pi-angle)*2;
s=linspace(0,log(tan(.25*coef)+sec(.25*coef))/coef,N);
se=linspace(0,log(tan((.25-eps)*coef)+sec((.25-eps)*coef))/coef,N);
x1 = 2*atan2(exp(coef*s)-1,exp(coef*s)+1)/coef;
xe = 2*atan2(exp(coef*se)-1,exp(coef*se)+1)/coef;

y1 = log(cos(coef*x1))*1/(coef)-T*coef;
