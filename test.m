


clear all;
load('data1')
sigma_d=max(cpg0x2Dsites1(:))/100;%domain sigma
sampling_d=sigma_d;

derived_sigma=sigma_d/sampling_d;

xi=round(cpg0x2Dsites1/sampling_d)+1;
max_x=max(xi);
numerator=zeros(1,max_x);
denominator=zeros(1,max(xi));


kernel_width=2*derived_sigma+1;

halfKernelWidth = floor( kernel_width / 2 );
kernel=[0 : kernel_width - 1];
kernel=kernel-halfKernelWidth;
kernel = (kernel.^2) / ( derived_sigma * derived_sigma ) ;
kernel = exp( -0.5 * kernel );%gaussian

for i=1:max_x
    %i
    mask=(xi==i);
    for col=1:size(methylation1,2)
        
        numerator(i,col)=sum(methylation1(mask,col));
        denominator(i,col)=sum(coverage1(mask,col));
        
    end
    
    
end
for col=1:size(methylation1,2)

    numerator(:,col)=conv(numerator(:,col),kernel,'same');
    denominator(:,col)=conv(denominator(:,col),kernel,'same');
    mask=(denominator(:,col)==0);
    numerator(mask,col)=0;
    denominator(mask)=1;
   
    smoothed_signal(:,col)=numerator(:,col)./denominator(:,col); 
    yi(:,col) = interp1([1:max_x],smoothed_signal(:,col),(cpg0x2Dsites1/sampling_d)+1);%%interpolation

end
%%coarse scale
figure;plot(smoothed_signal)
%%interpolated result
figure;plot(yi)