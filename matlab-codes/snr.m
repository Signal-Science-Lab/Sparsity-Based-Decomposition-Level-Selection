% Script to add noise to model data w/ SNR = 10, 50, 100
% ------------------------------------------------------ 
while abs(SNR - 10) > 0.0001
    sig = zeros(1,1024);
    noisevec = awgn(sig, 31);
    model = x + noisevec;
    SNR = max(x)/rms(model(1:200));
end
new_model_10 = model; % model + noise: SNR = 10
disp(fprintf('SNR = %4.3f', max(x)/rms(new_model_10(1:200)))) % verify SNR = 10

while abs(SNR - 50) > 0.0001
    sig = zeros(1,1024);
    noisevec = awgn(sig, 45);
    model = x + noisevec;
    SNR = max(x)/rms(model(1:200));
end
new_model_50 = model; % model + noise: SNR = 50
disp(fprintf('SNR = %4.3f', max(x)/rms(new_model_50(1:200)))) % verify SNR = 50

while abs(SNR - 100) > 0.0001
    sig = zeros(1,1024);
    noisevec = awgn(sig, 51);
    model = x + noisevec;
    SNR = max(x)/rms(model(1:200));
end
new_model_100 = model; % model+noise: SNR = 100
disp(fprintf('SNR = %4.3f', max(x)/rms(new_model_100(1:200)))) % verify SNR = 100
