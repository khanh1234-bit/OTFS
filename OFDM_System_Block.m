function [bit_errors] = OFDM_System_Block(N_subcarriers, N_syms, M_mod, velocity_kmh, SNR_val, fc, fs)
    % Tham số OFDM
    CP_len = N_subcarriers / 4; % Cyclic Prefix (thường là 1/4 hoặc 1/8)
    M_bits = log2(M_mod);
    
    % 1. Tạo dữ liệu
    num_bits = N_subcarriers * N_syms * M_bits;
    data_bits = randi([0, 1], num_bits, 1);
    data_syms = qammod(bi2de(reshape(data_bits, [], M_bits)), M_mod, 'gray');
    x_freq = reshape(data_syms, N_subcarriers, N_syms);
    
    % 2. Điều chế OFDM (IFFT + CP)
    x_time = ifft(x_freq, N_subcarriers) * sqrt(N_subcarriers);
    x_time_cp = [x_time(N_subcarriers-CP_len+1:end, :); x_time];
    s_tx = x_time_cp(:); % Chuỗi phát nối tiếp
    
    % 3. Kênh truyền (LTV Channel - Linear Time Variant)
    % Tái sử dụng logic tạo kênh để lấy thông số vật lý
    [taps, delay_taps, ~, chan_coef, max_doppler_Hz] = Channel_Gen_Velocity(1, 1, velocity_kmh, fc, fs);
    
    % Áp dụng kênh thay đổi theo thời gian
    % r(t) = sum( h_i * s(t - tau_i) * exp(j*2*pi*nu_i*t) )
    noise_var = 1 / SNR_val;
    r_rx = zeros(size(s_tx));
    
    % Tạo Doppler shifts cho từng tap (Hz)
    doppler_shifts = max_doppler_Hz * cos(rand(1, taps) * 2 * pi); 
    
    % Chiều dài tín hiệu
    L_sig = length(s_tx);
    time_idx = (0:L_sig-1).'/fs;
    
    for i = 1:taps
        % Tín hiệu trễ
        s_delayed = [zeros(delay_taps(i), 1); s_tx(1:end-delay_taps(i))];
        % Thêm Doppler pha quay
        path_rx = s_delayed .* chan_coef(i) .* exp(1j * 2 * pi * doppler_shifts(i) * time_idx);
        r_rx = r_rx + path_rx;
    end
    
    % Thêm nhiễu AWGN
    noise = sqrt(noise_var/2) * (randn(size(r_rx)) + 1j*randn(size(r_rx)));
    r_rx = r_rx + noise;
    
    % 4. Giải điều chế OFDM
    r_mat = reshape(r_rx, N_subcarriers + CP_len, N_syms);
    r_no_cp = r_mat(CP_len+1:end, :); % Bỏ CP
    y_freq = fft(r_no_cp, N_subcarriers) / sqrt(N_subcarriers);
    
    % 5. Equalizer (Zero Forcing đơn giản)
    % OFDM giả định kênh phẳng trong 1 symbol. Ta tính H tại tần số f
    % H_freq = FFT(h_time_impulses)
    % Lưu ý: Ở tốc độ cao, H thay đổi TRONG symbol, nhưng Rx chuẩn không biết điều đó.
    % Ta dùng đáp ứng kênh trung bình hoặc tại thời điểm đầu (để mô phỏng lỗi thực tế)
    
    H_est = zeros(N_subcarriers, 1);
    for i = 1:taps
        H_est = H_est + chan_coef(i) * exp(-1j * 2 * pi * (0:N_subcarriers-1)' * delay_taps(i) / N_subcarriers);
    end
    
    % Cân bằng kênh (Chia cho H_est)
    x_est = y_freq ./ repmat(H_est, 1, N_syms);
    
    % 6. Demapping & Tính lỗi
    data_demod = qamdemod(x_est(:), M_mod, 'gray');
    bits_est = reshape(de2bi(data_demod, M_bits), [], 1);
    bit_errors = sum(xor(bits_est, data_bits));
end