function xlinc=getlinearsegment(xxtr,noisex,n)
% CHEck Initial and Final Samples
   if  and( xxtr(1)>std(noisex) , xxtr(end)>-std(noisex))
        % 1st and Final Sample Test above noise
        if n==1
% %                     if length(xxtr)>2 % straight line
% %                         mslope=(xxtr(end)-xxtr(1))/length(xxtr);
% %                         xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
% %                         Apeaks=findpeaks(xxtr-xlinc);
% %                         Avalle=findpeaks(xlinc-xxtr);
% %                     else
% %                         Apeaks=[];
% %                         Avalle=[];
% %                     end
            %if isempty(Apeaks) % exp FIT
            if max(xxtr)==xxtr(1) % exp FIT
                xexp=fit([1:length(xxtr)]',xxtr','exp1');
                xlinc=xexp(1:length(xxtr))';
            else                % Zero Linea
                mslope=0;
                xlinc=mslope*([1:length(xxtr)]-1);
            end
            % Exponential  Decaying at the Start
        else
            mslope=0;
            xlinc=mslope*([1:length(xxtr)]-1);
            % Possible Calcium Transient
        end
    elseif and(xxtr(1)<-std(noisex),xxtr(end)<-std(noisex))
        % 1st Sample below Noise & Final Sample below noise
        %if and( numel(xxtr(xxtr>max([xxtr(1),xxtr(end)])))>numel(xxtr(xxtr<max([xxtr(1),xxtr(end)]))),...
        %    or(xxtr(1)>std(noisex),xxtr(end)>-std(noisex)) )
        if numel(xxtr(xxtr>max([xxtr(1),xxtr(end)])))>numel(xxtr(xxtr<max([xxtr(1),xxtr(end)])))
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            if  mslope>0
                mslope=0;
            end
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            disp('Ca2+ Transient')
        else
            xlinc=xxtr;
        end
    elseif and(xxtr(1)<-std(noisex),xxtr(end)>-std(noisex))
        % 1st Sample below & Final Sample above noise
        if or( numel(xxtr(xxtr>std(noisex)))>numel(xxtr(xxtr<0)),...
            or(xxtr(1)>std(noisex),xxtr(end)>-std(noisex)) )
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            disp('Ca2+ Transient')
        else
            xlinc=xxtr;
        end
    else
        if xxtr(end)>std(noisex)
            mslope=0;
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            disp('Ca2+ Transient')
        else
            % Otherwise Line
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            if sum(xlinc>xxtr)>numel(xxtr)/2
                % If line is above data, use smooth:
                xlinc=xxtr;
            end
        end
    end
end