% 3D Audio 
% Crosstalk Cancellation 
% Made by Yonghee Park & Zhiguang Eric Zhang

function individualized_crosstalk_cancellation(signame, ir30L, ir30R, ir330L, ir330R, filename)

    sig = audioread(signame);
    IR_30L = audioread(ir30L);
    IR_30R = audioread(ir30R);
    IR_330L = audioread(ir330L);
    IR_330R = audioread(ir330R);
    
    sig_len = length(sig);
    IR_len = length(IR_30L);
    total_len = sig_len + IR_len - 1;
    
    zero_sig = zeros(total_len - sig_len, 1);
    zero_IR = zeros(total_len - IR_len, 1);
    
    %stereo -> L and R 
    sig_L = [sig(:,1); zero_sig];
    sig_R = [sig(:,2); zero_sig];
    
    H30L = [IR_30L; zero_IR];
    H30R = [IR_30R; zero_IR];
    H330L = [IR_330L; zero_IR];
    H330R = [IR_330R; zero_IR];
    
    %fft
    sig_L = fft(sig_L);
    sig_R = fft(sig_R);
    
    H30L = fft(H30L);
    H30R = fft(H30R);
    H330L = fft(H330L);
    H330R = fft(H330R);
    
    %Crosstalk cancellation signal
    crosstalk_IR_L = (-1) .* (H330R./H330L);
    crosstalk_IR_R = (-1) .* (H30L./H30R);
    
    cross_cancel_L = sig_R .* crosstalk_IR_R;
    cross_cancel_R = sig_L .* crosstalk_IR_L;
    
    left_ch = sig_L + cross_cancel_L;
    right_ch = sig_R + cross_cancel_R;
    
    proc_IR_L = 1./(1-(H330R./H330L).^2);
    proc_IR_R = 1./(1-(H30L./H30R).^2);
    
    left_ch = left_ch .* proc_IR_L;
    right_ch = right_ch .* proc_IR_R;
    
    eq_IR_L = 1./H330L;
    eq_IR_R = 1./H30R;
    
    left_ch = left_ch .* eq_IR_L;
    right_ch = right_ch .* eq_IR_R;
    
    left_ch = ifft(left_ch);
    right_ch = ifft(right_ch);
    
    %Normalized
    left_ch = 0.999 * left_ch / max(abs(left_ch));
    right_ch = 0.999 * right_ch / max(abs(right_ch));
    
    sig_final = [left_ch, right_ch];
    sig_final = sig_final(1:sig_len,:);
    
    %Save
    wavwrite(sig_final, 44100, filename);

end