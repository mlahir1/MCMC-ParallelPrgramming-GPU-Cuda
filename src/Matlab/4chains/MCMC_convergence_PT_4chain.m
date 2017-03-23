% This m file takes 1 minute runtime on a Core i5 8 GB RAM machine
clc
clear
% Five iid Guassian Functions
f1=@(x1) 0.5*(exp(-(((x1-4)/1).^2)/2))/1/sqrt(2*pi)+ 0.5*(exp(-(((x1-6)/1).^2)/2))/1/sqrt(2*pi);
f2=@(x2) 0.5*(exp(-(((x2-5)/1).^2)/2))/1/sqrt(2*pi)+ 0.5*(exp(-(((x2-7)/1).^2)/2))/1/sqrt(2*pi);


f=@(x1,x2) (f1(x1).*f2(x2));

 expectancy= @(x1,x2) (x1*x2);
 
% MCMC method
accepted=0;
rejected=0;

accepted2=0;
rejected2=0;

accepted3=0;
rejected3=0;


accepted4=0;
rejected4=0;

exchanged=0;
EX_MCMC = 0;

% initialization
y1_now=1;y2_now=1;
y1_now2=1;y2_now2=1;
y1_now3=1;y2_now3=1;
y1_now4=1;y2_now4=1;


i=1;

VAR_TUNE1=0.458;
VAR_TUNE2=0.458;
VAR_TUNE3=0.458;
VAR_TUNE4=0.458;


while (1)
y1_c  = normrnd (y1_now, VAR_TUNE1);
y2_c  = normrnd (y2_now, VAR_TUNE1);

y1_c2  = normrnd (y1_now2, VAR_TUNE2);
y2_c2  = normrnd (y2_now2, VAR_TUNE2);

y1_c3  = normrnd (y1_now3, VAR_TUNE3);
y2_c3  = normrnd (y2_now3, VAR_TUNE3);

y1_c4  = normrnd (y1_now4, VAR_TUNE4);
y2_c4  = normrnd (y2_now4, VAR_TUNE4);

a1 = f(y1_c,y2_c)/f(y1_now,y2_now);
if rand < a1
    y1_now = y1_c;
    y2_now = y2_c;

    accepted=accepted+1;
else  
    rejected=rejected+1;
end


    
a2 = ((f(y1_c2,y2_c2)).^(9/16)) /((f(y1_now2,y2_now2)).^(9/16));
if rand < a2

	y1_now2 = y1_c2;
    y2_now2 = y2_c2;
    
    accepted2=accepted2+1;
else 
    rejected2=rejected2+1;
end


a3 = ((f(y1_c3,y2_c3)).^(4/16)) /((f(y1_now3,y2_now3)).^(4/16));
if rand < a3

	y1_now3 = y1_c3;
    y2_now3 = y2_c3;
    
    accepted3=accepted3+1;
else 
    rejected3=rejected3+1;
end

    
a4 = ((f(y1_c4,y2_c4)).^(1/16)) /((f(y1_now4,y2_now4)).^(1/16));
if rand < a4

	y1_now4 = y1_c4;
    y2_now4 = y2_c4;
    
    accepted4=accepted4+1;
else 
    rejected4=rejected4+1;
end



if accepted>rejected   VAR_TUNE1 = VAR_TUNE1*1.1; else VAR_TUNE1=VAR_TUNE1*0.9; end
if accepted2>rejected2 VAR_TUNE2 = VAR_TUNE2*1.1; else VAR_TUNE2=VAR_TUNE2*0.9; end
if accepted3>rejected3 VAR_TUNE3 = VAR_TUNE3*1.1; else VAR_TUNE3=VAR_TUNE3*0.9; end
if accepted4>rejected4 VAR_TUNE4 = VAR_TUNE4*1.1; else VAR_TUNE4=VAR_TUNE4*0.9; end


p1 = f(y1_now,y2_now).^(16/16);
p2 = f(y1_now2,y2_now2).^(9/16);
p1_2 = f(y1_now2,y2_now2).^(16/16);
p2_2 = f(y1_now,y2_now).^(9/16);
   if rand < (p2_2*p1_2)/(p1*p2)  
       exchange1 = y1_now;
       exchange2 = y2_now;

       y1_now  = y1_now2;
       y2_now  = y2_now2;

       y1_now2 = exchange1;
       y2_now2 = exchange2;
% 	   exchanged=exchanged+1;
   end
   

p3 = f(y1_now3,y2_now3).^(4/16);
p4 = f(y1_now4,y2_now4).^(1/16);
p3_2 = f(y1_now4,y2_now4).^(4/16);
p4_2 = f(y1_now3,y2_now3).^(1/16);
   if rand < (p4_2*p3_2)/(p3*p4)  
       exchange1 = y1_now3;
       exchange2 = y2_now3;

       y1_now3  = y1_now4;
       y2_now3  = y2_now4;

       y1_now4 = exchange1;
       y2_now4 = exchange2;
% 	   exchanged=exchanged+1;
   end


p3 = f(y1_now3,y2_now3).^(4/16);
p1 = f(y1_now,y2_now).^(16/16);
p3_2 = f(y1_now,y2_now).^(4/16);
p1_2 = f(y1_now3,y2_now3).^(16/16);
   if rand < (p1_2*p3_2)/(p1*p3)  
       exchange1 = y1_now3;
       exchange2 = y2_now3;

       y1_now3  = y1_now;
       y2_now3  = y2_now;

       y1_now = exchange1;
       y2_now = exchange2;
% 	   exchanged=exchanged+1;
   end

   
   

    EX_MCMC=EX_MCMC+expectancy(y1_now,y2_now);
    i=i+1;
       
    if ~mod (i,1000)
    EX_MCMC/i
    accepted
    rejected
    accepted2
    rejected2
	
    end

end

