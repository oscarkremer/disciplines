
clear all;
close all;

N = 128; % lenght of the PRBS sequence
M = 8; % number of samples per bit

if M<1 
    fprintf('\t\t ERROR\t\t\n');
    fprintf('Number of samples per bit is less than 1\n');
    return;
end
for k=1:N
    if (uk(k) == -1)
        for i=1:M
        waveformTemp(k,i) = [-1];
        end
        waveform{k} = waveformTemp(k,:);
    elseif (uk(k) == 1)
        for i=1:M
        waveformTemp(k,i) = [1];
        end
        waveform{k} = waveformTemp(k,:);
    end  
end
PRBSwave = cell2mat(waveform);
h = plot(PRBSwave,'k');
set (h, 'LineWidth', 2); % for a width of n
axis([0,length(PRBSwave),-1.5,1.5]);
title('PRBS waveform');
xlabel('samples')
ylabel('PRBS +1 / -1 values')
% end of script