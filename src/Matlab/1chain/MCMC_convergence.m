% This m file takes 1 minute runtime on a Core i5 8 GB RAM machine
clc
clear
% Five iid Guassian Functions
f1=@(x1) (exp(-(((x1-1)/1).^2)/2))/1/sqrt(2*pi);
f2=@(x2) (exp(-(((x2-1)/1).^2)/2))/1/sqrt(2*pi);
f3=@(x3) (exp(-(((x3-1)/1).^2)/2))/1/sqrt(2*pi);
f4=@(x4) (exp(-(((x4-1)/1).^2)/2))/1/sqrt(2*pi);
f5=@(x5) (exp(-(((x5-1)/1).^2)/2))/1/sqrt(2*pi);
f6=@(x6) (exp(-(((x6-1)/1).^2)/2))/1/sqrt(2*pi);
f7=@(x7) (exp(-(((x7-1)/1).^2)/2))/1/sqrt(2*pi);
f8=@(x8) (exp(-(((x8-1)/1).^2)/2))/1/sqrt(2*pi);
f9=@(x9) (exp(-(((x9-1)/1).^2)/2))/1/sqrt(2*pi);
% f10=@(x10) (exp(-(((x10-1)/1).^2)/2))/1/sqrt(2*pi);
% f11=@(x11) (exp(-(((x11-1)/1).^2)/2))/1/sqrt(2*pi);
% f12=@(x12) (exp(-(((x12-1)/1).^2)/2))/1/sqrt(2*pi);
% f13=@(x13) (exp(-(((x13-1)/1).^2)/2))/1/sqrt(2*pi);
% f14=@(x14) (exp(-(((x14-1)/1).^2)/2))/1/sqrt(2*pi);
% f15=@(x15) (exp(-(((x15-1)/1).^2)/2))/1/sqrt(2*pi);

% f1=@(x1) 1;
% f2=@(x2) 1;
% f3=@(x3) 1;
% f4=@(x4) 1;
% f5=@(x5) 1;
% f6=@(x6) 1;
% f7=@(x7) 1;
% f8=@(x8) 1;
% f9=@(x9) 1;
f10=@(x10) 1;
f11=@(x11) 1;
f12=@(x12) 1;
f13=@(x13) 1;
f14=@(x14) 1;
f15=@(x15) 1;

f=@(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15) (f1(x1).*f2(x2).*f3(x3).*f4(x4).*f5(x5).*f6(x6).*f7(x7).*f8(x8).*f9(x9).*f10(x10).*f11(x11).*f12(x12).*f13(x13).*f14(x14).*f15(x15));

% expectancy= @(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15) (x1*x2*x3*x4*x5*x6*x7*x8*x9*x10*x11*x12*x13*x14*x15);
expectancy=   @(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15) (x1*x2*x3*x4*x5*x6*x7*x8*x9);

% MCMC method
accepted=0;
rejected=0;

EX_MCMC = 0;

% initialization
y1(1)=1;y2(1)=1;y3(1)=1;y4(1)=1;y5(1)=1;y6(1)=1;y7(1)=1;y8(1)=1;y9(1)=1;y10(1)=1;y11(1)=1;y12(1)=1;y13(1)=1;y14(1)=1;y15(1)=1;
tic

i=1;
criteria=0;
while ~(criteria>100000)
y1_c  = normrnd (y1(i), 0.458);
y2_c  = normrnd (y2(i), 0.458);
y3_c  = normrnd (y3(i), 0.458);
y4_c  = normrnd (y4(i), 0.458);
y5_c  = normrnd (y5(i), 0.458);
y6_c  = normrnd (y6(i), 0.458);
y7_c  = normrnd (y7(i), 0.458);
y8_c  = normrnd (y8(i), 0.458);
y9_c  = normrnd (y9(i), 0.458);
y10_c = normrnd (y10(i), 0.458);
y11_c = normrnd (y11(i), 0.458);
y12_c = normrnd (y12(i), 0.458);
y13_c = normrnd (y13(i), 0.458);
y14_c = normrnd (y14(i), 0.458);
y15_c = normrnd (y15(i), 0.458);

if rand < f(y1_c,y2_c,y3_c,y4_c,y5_c,y6_c,y7_c,y8_c,y9_c,y10_c,y11_c,y12_c,y13_c,y14_c,y15_c)/f(y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i),y8(i),y9(i),y10(i),y11(i),y12(i),y13(i),y14(i),y15(i))
    y1(i+1) = y1_c;
    y2(i+1) = y2_c;
    y3(i+1) = y3_c;
    y4(i+1) = y4_c;
    y5(i+1) = y5_c;
    y6(i+1) = y6_c;
    y7(i+1) = y7_c;
    y8(i+1) = y8_c;
    y9(i+1) = y9_c;
    y10(i+1) = y10_c;
    y11(i+1) = y11_c;
    y12(i+1) = y12_c;
    y13(i+1) = y13_c;
    y14(i+1) = y14_c;
    y15(i+1) = y15_c;
    
    accepted=accepted+1;
else
    y1(i+1) = y1(i);
    y2(i+1) = y2(i);  
    y3(i+1) = y3(i);  
    y4(i+1) = y4(i);  
    y5(i+1) = y5(i);  
    y6(i+1) = y6(i);  
    y7(i+1) = y7(i);  
    y8(i+1) = y8(i);  
    y9(i+1) = y9(i);  
    y10(i+1) = y10(i);  
    y11(i+1) = y11(i);  
    y12(i+1) = y12(i);  
    y13(i+1) = y13(i);  
    y14(i+1) = y14(i);  
    y15(i+1) = y15(i);  
    
    rejected=rejected+1;
end
    EX_MCMC=EX_MCMC+expectancy(y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i),y8(i),y9(i),y10(i),y11(i),y12(i),y13(i),y14(i),y15(i));
    i=i+1;
    
    if (EX_MCMC/i)<1.05 & (EX_MCMC/i)>0.95
        criteria=criteria+1;
    else
        criteria=0;
    end
    
    if mod (i,1000)
    EX_MCMC/i
    accepted
    rejected
    end

end


toc