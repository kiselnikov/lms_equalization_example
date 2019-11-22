%--------------------------------------------------------------------------
% Company: ETMC Exponenta                                                        
% Engineer: Kiselnikov Andrei                                                    
%                                                                       
% Revision 0.01 - File Created 21.11.2019                                        
% This example provided strictly for educational purposes!
%
% Brief : LMS equalizer script
%--------------------------------------------------------------------------

clear all
close all

% Parameters --------------------------------------------------------------
%Number of bits
N =  128;

% The data sequence formation
data = round(rand(N,1));

% Raised cosine filter
roll_off_factor = 0.2;
span_in_symbols = 10;
symbols_per_step = 8;

% Equalizer length
eq_len = 21;

%Convergence multiplier
mu = 0.1;
%--------------------------------------------------------------------------

% Rcos fir design ---------------------------------------------------------
nyq_filter =  comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          roll_off_factor, ...
  'FilterSpanInSymbols',    span_in_symbols, ...
  'OutputSamplesPerSymbol', symbols_per_step);
%--------------------------------------------------------------------------

% Channel symbols design --------------------------------------------------
reference_signal =  upfirdn(data, nyq_filter.coeffs.Numerator, symbols_per_step);
%--------------------------------------------------------------------------

% Channel design ----------------------------------------------------------
channel_grp_delay = 10;
channel_coeffs = [-0.00448802244256685, ...
                   0.00645719073272638, ...
                   0.0194125083505277,  ...
                   0.0337393246952826,  ...
                   0.0486610611188258,  ...
                   0.0633171304003086,  ...
                   0.0768253273648968,  ...
                   0.0883475673098707,  ...
                   0.0971535155864083,  ...
                   0.102676757570639,   ...
                   0.104558743492511,   ...
                   0.102676757570639,   ...
                   0.0971535155864083,  ... 
                   0.0883475673098707,  ...
                   0.0768253273648968,  ...
                   0.0633171304003086,  ...
                   0.0486610611188258,  ...
                   0.0337393246952826,  ...
                   0.0194125083505277,  ...
                   0.00645719073272638, ...
                   -0.00448802244256685];
channel_output = filter (channel_coeffs,1,reference_signal);
channel_output = channel_output(channel_grp_delay+1:end);
%--------------------------------------------------------------------------

% LMS equalization performing ---------------------------------------------
%Adaptive filter coefficient initial array
w = zeros (eq_len,1);

%Error behavior
err = zeros (length(reference_signal)-2*eq_len-channel_grp_delay,1);

% Equalized output
equalized_output = zeros(length(reference_signal)-eq_len-channel_grp_delay,1);

for i = eq_len : (length(reference_signal)- channel_output - eq_len) 
	u = channel_output(i:-1:i+1-eq_len);
    equalized_output(i+1-eq_len)= w' * u; 
    err(i+1-eq_len) = reference_signal(i) - equalized_output(i+1-eq_len);
	w = w + mu * u  *  err(i+1-eq_len);
end
%--------------------------------------------------------------------------

% Visualization -----------------------------------------------------------
% Channel impulse response vs Equalizer response
subplot (2,2,1);
plot(channel_coeffs, 'ko')
hold on
plot(w, 'r*')
legend('Channel response','Eqlizer response')
grid on
title("Estimated weights by "+N*symbols_per_step+"samples") ;
% Channel signal vs Equalized signal
subplot (2,2,2);
plot(channel_output(eq_len:end));
hold on
plot(equalized_output);
grid on
legend('Channel signal','Eqlized signal')
title("Channel signal vs Equalized signal");

%Error between channel signal and reference signal
subplot (2,2,3);
plot(movmean(abs(err),span_in_symbols));
grid on
title("Average error behavior");
%Error between channel signal and reference signal
subplot (2,2,4);
plot(movmean(abs(reference_signal(1:length(reference_signal)-channel_grp_delay)-channel_output),span_in_symbols));
grid on
title("Average error on the Equalizer input");
%--------------------------------------------------------------------------
