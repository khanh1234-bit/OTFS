clc; clear all; close all;

%% --- Simulation Parameters ---
% Thông số hệ thống cơ bản
fc = 4e9;       % 4 GHz Carrier
fs = 15e3 * 32; % Bandwidth (giả định subcarrier 15kHz, 32 subcarriers)

% Kích thước lưới OTFS (N delay x M doppler)
N = 32; M = 32;
M_mod = 4; % 4-QAM (QPSK)
M_bits = log2(M_mod);

% Danh sách vận tốc cần mô phỏng (km/h)
velocities = [30, 120, 500];

% Dải SNR (dB)
SNR_dB = 5:5:25; % Chạy từ 5 đến 25 dB
SNR_lin = 10.^(SNR_dB/10);

% Số lượng khung hình chạy Monte Carlo (Giảm số này nếu máy chạy chậm)
N_fram = 50; 

% Biến lưu kết quả
ber_otfs = zeros(length(velocities), length(SNR_dB));
ber_ofdm = zeros(length(velocities), length(SNR_dB));

fprintf('Bắt đầu mô phỏng so sánh OTFS vs OFDM...\n');

%% --- Main Loop ---
for v_idx = 1:length(velocities)
    vel = velocities(v_idx);
    fprintf('Running Velocity: %d km/h\n', vel);
    
    for s_idx = 1:length(SNR_dB)
        snr_val = SNR_lin(s_idx);
        sigma_2 = 1/snr_val; % Noise variance (assuming signal power norm to 1)
        
        total_err_otfs = 0;
        total_err_ofdm = 0;
        total_bits = 0;
        
        parfor ifram = 1:N_fram % Dùng parfor để chạy song song nếu có thể
            % --- 1. OTFS Simulation ---
            % Tạo bits
            bits_per_frame = N*M*M_bits;
            data_bits = randi([0,1], bits_per_frame, 1);
            data_syms = qammod(bi2de(reshape(data_bits, [], M_bits)), M_mod, 'gray');
            x = reshape(data_syms, N, M);
            
            % OTFS Mod
            s = OTFS_modulation(N, M, x);
            
            % Channel Gen (Velocity dependent)
            [taps, delay_taps, Doppler_taps, chan_coef, ~] = Channel_Gen_Velocity(N, M, vel, fc, fs);
            
            % OTFS Channel Output (Sử dụng hàm cũ của bạn)
            r = OTFS_channel_output(N, M, taps, delay_taps, Doppler_taps, chan_coef, sigma_2, s);
            
            % OTFS Demod
            y = OTFS_demodulation(N, M, r);
            
            % OTFS MP Detector (Sử dụng hàm cũ của bạn)
            % Lưu ý: MP Detector cần đúng tham số Doppler_taps (dạng số thực hoặc index)
            % Để code MP cũ hoạt động tốt với Doppler thực (lẻ), ta cần đảm bảo
            % logic trong MP detector xử lý được fractional Doppler.
            % Code cũ của bạn xử lý Doppler shift integer. 
            % Để đơn giản hóa cho demo này, ta làm tròn Doppler index trong MP detector
            % hoặc chấp nhận sai số nhỏ. Ở đây ta truyền Doppler_taps trực tiếp.
            x_est_syms = OTFS_mp_detector(N, M, M_mod, taps, delay_taps, Doppler_taps, chan_coef, sigma_2, y);
            
            % Calculate Errors
            data_demapping = qamdemod(x_est_syms, M_mod, 'gray');
            bits_est = reshape(de2bi(data_demapping, M_bits), bits_per_frame, 1);
            err_otfs = sum(xor(bits_est, data_bits));
            
            % --- 2. OFDM Simulation ---
            % Để công bằng, OFDM dùng cùng số lượng symbol và subcarrier (M subcarr, N time slots)
            err_ofdm = OFDM_System_Block(M, N, M_mod, vel, snr_val, fc, fs);
            
            % Accumulate
            total_err_otfs = total_err_otfs + err_otfs;
            total_err_ofdm = total_err_ofdm + err_ofdm;
            
        end
        
        % Tính BER trung bình
        ber_otfs(v_idx, s_idx) = total_err_otfs / (N_fram * N * M * M_bits);
        ber_ofdm(v_idx, s_idx) = total_err_ofdm / (N_fram * N * M * M_bits);
        
        fprintf('  SNR %d dB -> BER OTFS: %.5f | BER OFDM: %.5f\n', SNR_dB(s_idx), ber_otfs(v_idx, s_idx), ber_ofdm(v_idx, s_idx));
    end
end

%% --- Plotting Results ---
figure('Color', 'w');
markers = {'-o', '-s', '-d'};
colors_otfs = {'r', 'b', 'k'}; % Red, Blue, Black for OTFS
colors_ofdm = {'r--', 'b--', 'k--'}; % Dashed for OFDM

semilogy(SNR_dB, ber_otfs(1,:), [colors_otfs{1} markers{1}], 'LineWidth', 2, 'DisplayName', 'OTFS 30 kmph'); hold on;
semilogy(SNR_dB, ber_otfs(2,:), [colors_otfs{2} markers{2}], 'LineWidth', 2, 'DisplayName', 'OTFS 120 kmph');
semilogy(SNR_dB, ber_otfs(3,:), [colors_otfs{3} markers{3}], 'LineWidth', 2, 'DisplayName', 'OTFS 500 kmph');

semilogy(SNR_dB, ber_ofdm(1,:), [colors_ofdm{1} markers{1}], 'LineWidth', 1.5, 'DisplayName', 'OFDM 30 kmph');
semilogy(SNR_dB, ber_ofdm(2,:), [colors_ofdm{2} markers{2}], 'LineWidth', 1.5, 'DisplayName', 'OFDM 120 kmph');
semilogy(SNR_dB, ber_ofdm(3,:), [colors_ofdm{3} markers{3}], 'LineWidth', 1.5, 'DisplayName', 'OFDM 500 kmph');

grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER Performance: OTFS vs OFDM (Ideal Pulses)');
legend('Location', 'southwest');
ylim([1e-5 1]);