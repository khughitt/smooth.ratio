%{
qt_cov=quantile(coverage,0.98,1);

for c=1:size(coverage,2)   
    tmp=coverage(:,c);
    tmp(tmp>qt_cov(c))=qt_cov(c);
    coverage(:,c)=tmp;    
end
%}
bins=9;

radius=350;

smoothed_percentage=zeros(size(coverage));

for c=1:size(coverage,2)   

    c_cum=cumsum(coverage(:,c),1);
    m_cum=cumsum(methylation(:,c),1);
    smoothed_percentage(1:radius,c)=m_cum(1+radius:2*radius)./c_cum(1+radius:2*radius);
    smoothed_percentage(end-radius+1:end,c)=(repmat(m_cum(end),[radius,1])-m_cum(end-2*radius+1:end-radius))./(repmat(c_cum(end),[radius,1])-c_cum(end-2*radius+1:end-radius));
    smoothed_percentage(radius+1:end-radius,c)=(m_cum(1+2*radius:end)-m_cum(1:end-2*radius))./(c_cum(1+2*radius:end)-c_cum(1:end-2*radius));
end




