clc;
clear all;
close all;

N=2000; % signal length
M=11; % weight number
w0 = zeros([M 1]); % inital weight
run_step=200; % train step
% mu = [0.15 0.075 0.01 0.005]; % desired signal's delay time 
%% problem 1
tau=7; % delay time
W=[2.9,3.1,3.3,3.5]; % Channel distortion coefficient
mu = 0.01; % step size
for idx=1:4
    h=zeros([1 3]); % channel distortion
    for n=1:3
        h(n)=(1+cos((2*pi)/W(idx)*(n-2)))/2;
    end
    
    
    e=zeros([run_step N]);
    % run 200 independent signal
    for i=1:run_step
        w=w0;
        d=2*randi([0,1],1,N)-1;
        d_delay=[zeros([1,tau]) d(1:N-tau)];
        u = conv(d,h);
        U = []; % matrix of u
        for j=1:M
            U = [U;u(1:N)];
            u=[0 u];
        end
        for n=1:N
            y = w' * U(:,n); % y=w*U
            e(i,n)=d_delay(n)-y; % error = d(n) - y(n) 
            w = w + mu*conj(e(i,n))*U(:,n);
        end
        
    end
    figure(1);
    subplot(2,2,idx);
    plot(mean(e.^2)); 
    axis([0 2000 0 1]);
    xlabel('n');
    ylabel('e(n)^2');
    title("W="+num2str(W(idx)));
end

%% problem 2
tau=[7,9,11,15]; % delay time
W=2.9; % Channel distortion coefficient
mu = 0.01; % step size
for idx=1:4
    h=zeros([1 3]); % channel distortion
    for n=1:3
        h(n)=(1+cos((2*pi)/W*(n-2)))/2;
    end
    
    
    e=zeros([run_step N]);
    % run 200 independent signal
    for i=1:run_step
        w=w0;
        d=2*randi([0,1],1,N)-1;
        d_delay=[zeros([1,tau(idx)]) d(1:N-tau(idx))];
        u = conv(d,h);
        U = []; % matrix of u
        for j=1:M
            U = [U;u(1:N)];
            u=[0 u];
        end
        for n=1:N
            y = w' * U(:,n); % y=w*U
            e(i,n)=d_delay(n)-y; % error = d(n) - y(n) 
            w = w + mu*conj(e(i,n))*U(:,n);
        end
        
    end
    figure(1);
    subplot(2,2,idx);
    plot(mean(e.^2)); 
    axis([0 2000 0 2]);
    xlabel('n');
    ylabel('e(n)^2');
    title("τ ="+num2str(tau(idx)));
end
%% problem 3
tau=7; % delay time
W=2.9; % Channel distortion coefficient
mu = [0.3,0.1,0.05,0.01]; % step size
for idx=1:4
    h=zeros([1 3]); % channel distortion
    for n=1:3
        h(n)=(1+cos((2*pi)/W*(n-2)))/2;
    end
    
    
    e=zeros([run_step N]);
    % run 200 independent signal
    for i=1:run_step
        w=w0;
        d=2*randi([0,1],1,N)-1;
        d_delay=[zeros([1,tau]) d(1:N-tau)];
        u = conv(d,h);
        U = []; % matrix of u
        for j=1:M
            U = [U;u(1:N)];
            u=[0 u];
        end
        for n=1:N
            y = w' * U(:,n); % y=w*U
            e(i,n)=d_delay(n)-y; % error = d(n) - y(n) 
            w = w + mu(idx)*conj(e(i,n))*U(:,n);
        end
        
    end
    figure(1);
    subplot(2,2,idx);
    plot(mean(e.^2)); 
    axis([0 2000 0 2]);
    xlabel('n');
    ylabel('e(n)^2');
    title("µ ="+num2str(mu(idx)));
end
%% problem 4 - Noise 

tau=7; % delay time
W=2.9; % Channel distortion coefficient
mu = 0.01; % step size
snr=[10,5,3,1];
for idx=1:4
    h=zeros([1 3]); % channel distortion
    for n=1:3
        h(n)=(1+cos((2*pi)/W*(n-2)))/2;
    end
    
    
    e=zeros([run_step N]);
    % run 200 independent signal
    for i=1:run_step
        % signal
        w=w0;
        d=2*randi([0,1],1,N)-1;
        d_delay=[zeros([1,tau]) d(1:N-tau)];
        u = conv(d,h);
        u=u(1:N);
        % noise 
        noise = randn([1 N]);
        noise =noise -mean(noise);
        signal_power=1/N*sum(u.*u);
        noise_variance = signal_power / (10^(snr(idx)/10));
        noise = (sqrt(noise_variance) / std(noise)) *noise;
        
        
        u = u+ noise;
        tmp=u;
        U = []; % matrix of u
        for j=1:M
            U = [U;u(1:N)];
            u=[0 u];
        end
        for n=1:N
            y = w' * U(:,n); % y=w*U
            e(i,n)=d_delay(n)-y; % error = d(n) - y(n) 
            w = w + mu*conj(e(i,n))*U(:,n);
        end
        
    
    end
    figure(1);
    subplot(2,2,idx);
    plot(mean(e.^2)); 
    axis([0 2000 0 2]);
    xlabel('n');
    ylabel('e(n)^2');
    title("white noise SNR :"+num2str(snr(idx)));
end

