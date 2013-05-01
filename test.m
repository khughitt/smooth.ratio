% Box Smooth Test
% Chengxi Ye <cxy@cs.umd.edu>
% 2013/04/29
clear all;
load('data/sample')

sigma_d    = (max(cpgsites(:)) - min(cpgsites(:))) / 100; % Spatial domain
sampling_d = sigma_d;

derived_sigma=sigma_d / sampling_d;

xi    = round((cpgsites - min(cpgsites)) / sampling_d) + 1;
max_x = max(xi);
numerator   = zeros(max_x, size(methylation, 2));
denominator = zeros(max_x, size(methylation, 2));

kernel_width = 2 * derived_sigma + 1;

halfKernelWidth = floor(kernel_width / 2);
kernel = [0:kernel_width - 1];
kernel = kernel - halfKernelWidth;
kernel = (kernel.^2) / (derived_sigma * derived_sigma);
kernel = exp(-0.5 * kernel); %gaussian

for i=1:max_x
    mask=(xi==i);
    for col=1:size(methylation,2)
        numerator(i, col)   = sum(methylation(mask, col));
        denominator(i, col) = sum(coverage(mask, col));
    end
end

for col=1:size(methylation,2)
    numerator(:,col)   = conv(numerator(:,col),kernel,'same');
    denominator(:,col) = conv(denominator(:,col),kernel,'same');
    mask = (denominator(:,col)==0);
    numerator(mask,col)   = 0;
    denominator(mask,col) = 1;
   
    smoothed_signal(:,col) = numerator(:,col) ./ denominator(:,col);
    
    % interpolation
    yi(:,col) = interp1([1:max_x], smoothed_signal(:,col), 
                        (cpgsites/sampling_d) + 1);
end

%%coarse scale
figure;plot(smoothed_signal)

%%interpolated result
%%figure;plot(yi)
