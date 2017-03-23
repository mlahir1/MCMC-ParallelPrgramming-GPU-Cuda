% This m file takes 1 minute runtime on a Core i5 8 GB RAM machine
clc
clear
% Five iid Guassian Functions
f1=@(x1) 0.5*(exp(-(((x1+18)/1).^2)/2))/1/sqrt(2*pi)+ 0.5*(exp(-(((x1-18)/1).^2)/2))/1/sqrt(2*pi);
f2=@(x2) 0.5*(exp(-(((x2+14)/1).^2)/2))/1/sqrt(2*pi)+ 0.5*(exp(-(((x2-14)/1).^2)/2))/1/sqrt(2*pi);


f=@(x1,x2) (f1(x1).*f2(x2));

 expectancy= @(x1,x2) (x1*x2);
 
% MCMC method
accepted=0;
rejected=0;

accepted2=0;
rejected2=0;

exchanged=0;
EX_MCMC = 0;

% initialization
y1_now=1;y2_now=1;
y1_now2=1;y2_now2=1;


i=1;

VAR_TUNE1=0.458;
VAR_TUNE2=0.458;

while (1)
y1_c  = normrnd (y1_now, VAR_TUNE1);
y2_c  = normrnd (y2_now, VAR_TUNE1);

y1_c2  = normrnd (y1_now2, VAR_TUNE2);
y2_c2  = normrnd (y2_now2, VAR_TUNE2);

a1 = f(y1_c,y2_c)/f(y1_now,y2_now);
if rand < a1
    y1_now = y1_c;
    y2_now = y2_c;

    accepted=accepted+1;
else  
    rejected=rejected+1;
end



    
a2 = ((f(y1_c2,y2_c2)).^(1/17)) /((f(y1_now2,y2_now2)).^(1/17));
if rand < a2

	y1_now2 = y1_c2;
    y2_now2 = y2_c2;
    
    accepted2=accepted2+1;
else 
    rejected2=rejected2+1;
end


if accepted>rejected   VAR_TUNE1 = VAR_TUNE1*1.1; else VAR_TUNE1=VAR_TUNE1*0.9; end
if accepted2>rejected2 VAR_TUNE2 = VAR_TUNE2*1.1; else VAR_TUNE2=VAR_TUNE2*0.9; end


p1 = f(y1_now,y2_now);
p2 = f(y1_now2,y2_now2);
p1_2 = p1.^(1/17);
p2_2 = p2.^(1/17);


   if rand < (p2*p1_2)/(p1*p2_2)  
       exchange1 = y1_now;
       exchange2 = y2_now;

       y1_now  = y1_now2;
       y2_now  = y2_now2;

       y1_now2 = exchange1;
       y2_now2 = exchange2;
	   exchanged=exchanged+1;
	   
   end


    EX_MCMC=EX_MCMC+expectancy(y1_now,y2_now);
    i=i+1;
       
    if ~mod (i,1000)
    EX_MCMC/i
    accepted
    rejected
    accepted2
    rejected2
	exchanged	
    end

end

