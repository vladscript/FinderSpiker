%% Clustering Coefficient
function [CoefClustIdx h]=GetClusteringCoef(Tree,Name,ActiveDataPeaks,AdjacentPeaks,numFig)
    %AdjacentPeaks=GetAdjacencyFromPeaks(handles,ActiveDataPeaks);
    %Method2Plot=get(handles.Clust2PlotPopupmenu,'Value');
    result=zeros(10);
    %disp('Matriz adyacente')
    %disp(AdjacentPeaks)
    for i=1:10,
        
        idxHC=cluster(Tree,Name,i);
        
        %switch Method2Plot
         %   case 1
                idxMethod=idxHC;
          %  case 2
           %     idxMethod=idxKmeans;
            %case 3
             %   idxMethod=idxHCRed;
           % case 4
            %    idxMethod=idxKmeansRed;
        %end
        %[CellsGroups1 CellsActGroups1]=GetCellsGroup(ActiveDataPeaks,idxMethod);
        
        Groups=max(idxMethod);
        Cells=size(ActiveDataPeaks,2);
        CellsGroups=zeros(Groups,Cells);
        CellsActGroups1=zeros(Groups,Cells);
        for ii=1:Groups
            PeaksGroupIdx= idxMethod==ii;
            CellsActGroups1(ii,:)=sum(ActiveDataPeaks(PeaksGroupIdx,:),1);
            CellsGroup= logical(CellsActGroups1(ii,:));
            CellsGroups(ii,CellsGroup)=ii+1;
        end
        
        %disp('Grupos')
        %disp(CellsActGroups1)
        [m,n]=size(CellsActGroups1);%m=numero de grupos y n=numero de elementos
        neigNode=zeros(1,n);
        neuGroup=zeros(n);
        for j=1:m, %j=grupo
            sumCoef=0;
            numElem=0;
            for k=1:n,
                 if(CellsActGroups1(j,k)>0)
                     neuGroup(k)=k;       
                 else
                     neuGroup(k)=0; 
                 end
            end
            newNeuGroup=find(neuGroup);
            for k=1:n,%k=numero de neurona
                sumLink=0;     
                sumElem=1;
                if(CellsActGroups1(j,k)>0)%la neurona k pertenece al grupo j
                    for l=1:n,
                        ind=find(newNeuGroup==l);
                        if(AdjacentPeaks(k,l)>0 & ind>0)
                        %if(AdjacentPeaks(k,l)>0)
                            neigNode(1,sumElem+1)=l;
                            sumElem=sumElem+1;
                        end
                    end
                    %calculo de la e
                    for x1=1:sumElem
                        for x2=1:sumElem
                            if(AdjacentPeaks(x1,x2)>0)
                                sumLink=sumLink+1;
                            end
                        end
                    end
                    %disp(' ');
                    %disp(sumLink);
                    %disp(sumElem*(sumElem-1));                    
                    sumCoef=sumCoef+sumLink/(sumElem*(sumElem-1));
                    numElem=numElem+1;
                    
                end
            end
            result(i,j)=sumCoef/numElem;
        end
        
    end
    meanGroups=zeros(1,10);
    for i=1:10,
        for j=1:i,
            meanGroups(i)=meanGroups(i)+result(i,j);
        end
        meanGroups(i)=meanGroups(i)/i;
    end
    for i=1:10,
        standD(i)=std(result(i,1:i));
    end
    for i=1:10,
        CoefClustIdx(i)= meanGroups(i)/standD(i);
    end
    CoefClustIdx(1)=0;
    %disp('Desviaciones Estandar');
    %disp(standD);    
    %disp('Media');
    %disp(meanGroups);
    %disp('Coeficiente');
    %disp(CoefClustIdx);
    h=max(CoefClustIdx);


    figure(numFig);
    plot(CoefClustIdx,'-b');
    hold on
    [hIdx h]=max(CoefClustIdx);
    plot(h,hIdx,'*r')
    title(['Clustering Coefficient (' num2str(h) ' groups recommended)'])
end