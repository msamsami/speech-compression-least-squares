clear
close all
clc

%STEP 1--------------------------------------------------
[y,Fs] = audioread('failure1FIX.wav');
sound(y,Fs);
pause(2.5);
figure('Name','y1 signal range','NumberTitle','off')
plot (y);
%STEP 1--------------------------------------------------


%STEP 2--------------------------------------------------
q = (max(y)-min(y))/10;   % L=10
yq = round(y/q)*q;
figure('Name','quantized signal range (L=10)','NumberTitle','off')
plot (yq);

sound(yq,Fs);
pause(2.5);

ly = length(y);  % Length of the original signal matrix (y)
for i = 1:16
    q=(max(y)-min(y))/i;
    yq=round(y/q)*q;
    MSE(1,i) = (y-yq)'*(y-yq)/ly;
end

figure('Name','MSE 1','NumberTitle','off')
plot (MSE);
%STEP 2--------------------------------------------------


%STEP 3 and 4--------------------------------------------

%################
% In this part we add some zeros (when necessary) at the beginning
% of the matrix y, in order to make its size divisible by 160 :
if mod(ly,160) ~= 0
    ydd = [zeros(mod(ly,160), 1); y];
end
%################

%################
% In this part we add 10 zeros at the first of the
% matrix y and make a new matrix called shY :
shY = [zeros(10,1); y];
%################

%################
% In this part we make the arguments of the toeplitz function :

y160 = y(1:160, 1);
y1 = shY(10:169, 1)';

for i = 1:10 
    y2(11-i) = shY(i,1);
end
%################

%################
% In this part we arrange the equation Aa=y and solve it to
% find a(k)s and also e(n)s :
AA = toeplitz(y1,y2); % A
aa = AA\y160;        % a(k) coefficients
error = y160-AA*aa;   % e(n) errors
%################

aa = [];
error = [];
for k = 1:ly/160
    y1=[];
    y2=[];
    err=[];
    bb=[];
    
    %################
    % Now we find a(k)s and e(n)s for all blocks of data : 
    y160 = y((k-1)*160+1:(k-1)*160+160, 1);
    y1 = shY((k-1)*160+10:(k-1)*160+169, 1)';

    for i = 1:10 
        y2(11-i) = shY((k-1)*160+i,1);
    end

    AA = toeplitz(y1,y2); 
    bb = AA\y160;
    err = y160-AA*bb;
    
    %################
    % In this part we add the a(k)s and e(n)s of this block to the
    % matrix related to all a(k)s and all e(n)s :
    aa = [aa; bb];
    error = [error; err];
    %################
end

figure('Name','Error 1')
plot (error);

%STEP 3 and 4--------------------------------------------



%STEP 5--------------------------------------------------
lyn = 0;  % initial length of the new reconstructed matrix (ynew)

for k = 0:ly/160 - 1
    y160 = 0;
    bb = 0;
    
    for i=1:10
        bb(i,1) = aa(k*10 + i,1);
    end
    
    for i=1:160
        y160(i,1) = error(k*160 + i, 1);
        for j=1:10
            y160(i,1) = y160(i,1) + bb(j,1) * shY(k*160 + 10 + i - j,1);
        end
    end
    
    for j=1:160
        ynew(lyn + j, 1) = y160(j,1);
    end
    
    lyn = length(ynew);
end

% ################
% We see that ynew is exactly the same as the original signal matrix y
% ################

sound(ynew, Fs);
pause(2.5);

%STEP 5--------------------------------------------------


%STEP 6--------------------------------------------------

q = (max(error)-min(error))/10;   %L=10
errorq = round(error/q)*q;

yhatr = zeros(10,1);
yhatyj = zeros(10,1);

lyh = 10; 

for k = 0:ly/160 - 1
    x=0;
    bb=[];
    
    bb = aa(k*10+1:k*10+10, 1);
    
    for i=1:160
        x = errorq(k*160 + i, 1);
        yj = 0;
        for j = 1:10
            x = x + bb(j,1) * yhatr(k*160 + 10 + i - j, 1);
            yj = yj + bb(j,1) * yhatr(k*160 + 10 + i - j, 1);
        end
        yhatr(lyh + i, 1) = x;
        yhatyj(lyh + i, 1) = yj;
    end
    lyh = length(yhatr);
end

yhat = yhatr(11:end, 1);
yhatfff = yhatyj(11:end, 1);


figure('Name','yhat1 signal','NumberTitle','off')
plot(yhat);
sound(yhat, Fs);
pause(2.5);

for d = 1:16
    yhatr = [];
    yhat = [];
    errorq = [];
    err = [];
    errq = [];
    leq = 0;

    for j = 0:(length(error)/160)-1
        err = error(j*160+1:j*160+160);

        q = (max(err)-min(err))/d;
        errq = round(err/q)*q;
        
        errorq = [errorq; errq(1:160, 1)];  
    end
    
    leq = length(errorq);
    yhatr = zeros(10,1);
    lyh=10; 
    
    for k = 0:ly/160 - 1
        x = 0;
        bb = [];

        bb = aa(k*10+1:k*10+10, 1);

        for i = 1:160
            x = errorq(k*160 + i, 1);
            for j = 1:10
                x= x + bb(j,1) * yhatr(k*160 + 10 + i - j, 1);
            end
            yhatr(k*160 + 10 + i, 1)=x;
        end
        
    end

    yhat = yhatr(11:end, 1);
    MSEhat(1,d)=(y-yhat)'*(y-yhat)/ly;
end

figure('Name','MSEhat 1','NumberTitle','off')
plot (MSEhat);

%STEP 6--------------------------------------------------
