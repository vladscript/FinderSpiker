function [x y xfit yfit slope R2] = FitPowerLaw(Adjacency)

    Links=sum(Adjacency,1);

    if ~Links
        return;
    end
    
    max_Links=max(Links);
    if ~max_Links
        return;
    end
    
    N=length(Links);
    nbins=round(log2(N)+1);
    [y x]=hist(Links,nbins);
    
    % Plot power law fit (only neurons with at least one link)
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    
end