function [significativeThModule]=testCoactivityGroup(raster,Ci,alpha)
    n=max(Ci);
    [N,F]=size(raster);
    significativeThModule=zeros(n,1);
    for i=1:n
        moduleNeurons=find(Ci==i);
        newC=length(moduleNeurons);
        newData=zeros(newC,F);
        for j=1:newC
            newData(j,:)=raster(moduleNeurons(j),:);
        end
        Co=sum(newData,1);
        disCumm=zeros(max(Co),1);
        for j=1:max(Co)
            disCumm(j)=length(Co(find(Co==j)));
        end
        sumCumm=0;
        sumTot=sum(disCumm);
        for j=1:max(Co)
            sumCumm=sumCumm+disCumm(j);
            disCumm(j)=sumCumm/sumTot;
        end
        num=find(disCumm>alpha,1,'first');
        if ~isempty(num)
            significativeThModule(i)=find(disCumm>alpha,1,'first');
        else
            significativeThModule(i)=0;
        end  
    end    
end