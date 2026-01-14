function [taps, delay_taps, Doppler_taps, chan_coef, max_doppler_Hz] = Channel_Gen_Velocity(N, M, velocity_kmh, fc, fs)
    %% 1. Tính toán vật lý
    c = 3e8; 
    % Tần số Doppler cực đại (Hz)
    max_doppler_Hz = (velocity_kmh / 3.6) * fc / c;
    
    % Độ phân giải Doppler của hệ thống (Hz per grid point)
    doppler_resolution = fs / (N * M); 
    
    % Quy đổi Doppler Hz sang số nguyên (Grid Index)
    % ceil để đảm bảo ít nhất là 1 nếu vận tốc > 0
    max_k = ceil(max_doppler_Hz / doppler_resolution);
    
    %% 2. Cấu hình Tap
    taps = 4; 
    
    % Delay taps (Cố định hoặc ngẫu nhiên số nguyên)
    delay_taps = [0, 1, 2, 3]; 
    
    % Doppler taps: Bắt buộc là SỐ NGUYÊN (Integer)
    % Chọn ngẫu nhiên trong khoảng [-max_k, max_k]
    if max_k == 0
        Doppler_taps = zeros(1, taps);
    else
        Doppler_taps = randi([-max_k, max_k], 1, taps);
    end
    
    % Hệ số kênh
    pow_prof = (1/taps) * ones(1, taps);
    chan_coef = sqrt(pow_prof) .* (sqrt(1/2) * (randn(1, taps) + 1i*randn(1, taps)));
end