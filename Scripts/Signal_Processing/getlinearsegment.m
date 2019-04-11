function xlinc=getlinearsegment(xxtr,StdNoise,n)
% CHEck Initial and Final Samples
mslope=(xxtr(end)-xxtr(1))/numel(xxtr);
xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
   if  and( xxtr(1)>StdNoise , xxtr(end)>-StdNoise)
        % 1st and Final Sample Test above noise
        if n==1
            %if isempty(Apeaks) % exp FIT
            if max(xxtr)==xxtr(1) || xxtr(1)>4*StdNoise  % exp FIT
                % Test if there is a peak in there:
                % 1st Derivative 
                ZeroCrosses=find(diff(sign(diff(xxtr))))+1;
                if numel(xxtr)>3
                    % Seek Valleys
                    [~,B]=findpeaks(-xxtr);
                else
                    B=1;
                end
                if isempty(B)
                    B=1;
                else
                    B=B(1);
                end
                % disp('CuRve3S5>>>')
                if ~isempty(ZeroCrosses)
                    mslope=(xxtr(end)-xxtr(B))/length(xxtr);
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(B);
                    xtest=xxtr-xlinc;
                    % what comes around->works-around
                    if numel(xtest(xtest<0))>numel(xtest(xtest>=0)) &&...
                            numel(xxtr<0)
                        % xexp=fit([1:length(xxtr)]',xxtr','exp1');
                        % xlinc=xexp(1:length(xxtr))';
                        xlinc=xxtr;
                    end
                else
                    xexp=fit([1:length(xxtr)]',xxtr','exp1');
                    xlinc=xexp(1:length(xxtr))';
                end
            else                % Zero Linea
                mslope=0;
                xlinc=mslope*([1:length(xxtr)]-1);
            end
            % Exponential  Decaying at the Start
        else
            % mslope=0;
            % xlinc=mslope*([1:length(xxtr)]-1);
            % Possible Calcium Transient
        end
    elseif and(xxtr(1)<0,xxtr(end)<-StdNoise)
    % 1st Sample below Noise & Final Sample below noise
        if numel(xxtr(xxtr>max([xxtr(1),xxtr(end)])))>numel(xxtr(xxtr<max([xxtr(1),xxtr(end)])))
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            % xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            xtest=mslope*([1:length(xxtr)]-1)+xxtr(1);
            if numel(xxtr)>2
                A=findpeaks(xxtr-xtest);
                if ~isempty(A(A>StdNoise))
                    xlinc=xtest;
                    disp('>> Ca2+ Possible Transient Rescued')
                else
                    if numel(xlinc)>3
                        PeakAmp=findpeaks(xxtr-xlinc);
                        PeakAmp=PeakAmp(PeakAmp>StdNoise);
                    else
                        PeakAmp=[];
                    end
                    if isempty(PeakAmp) && sum(xlinc>xxtr)>numel(xxtr)/2
                        % If line is above data, use smooth:
                        xlinc=xxtr;
                    end
                    disp('>>> Detrending >>>')
                end
            else
                if numel(xlinc)>3
                    PeakAmp=findpeaks(xxtr-xlinc);
                    PeakAmp=PeakAmp(PeakAmp>StdNoise);
                else
                    PeakAmp=[];
                end
                if isempty(PeakAmp) && sum(xlinc>xxtr)>numel(xxtr)/2
                    % If line is above data, use smooth:
                    xlinc=xxtr;
                end
                disp('WTF');
            end
            
        else
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xtest=mslope*([1:length(xxtr)]-1)+xxtr(1);
            if numel(xxtr)>2
                A=findpeaks(xxtr-xtest);
                if ~isempty(A(A>StdNoise))
                    xlinc=xtest;
                else
                    if numel(xtest)>3
                        PeakAmp=findpeaks(xxtr-xtest);
                        PeakAmp=PeakAmp(PeakAmp>StdNoise);
                    else
                        PeakAmp=[];
                    end
                    if isempty(PeakAmp) && sum(xtest>xxtr)>numel(xxtr)/2
                        % If line is above data, use smooth:
                        xlinc=xxtr;
                    end
                end
            else
                if numel(xlinc)>3
                    PeakAmp=findpeaks(xxtr-xlinc);
                    PeakAmp=PeakAmp(PeakAmp>StdNoise);
                else
                    PeakAmp=[];
                end
                if isempty(PeakAmp) && sum(xlinc>xxtr)>numel(xxtr)/2
                    % If line is above data, use smooth:
                    xlinc=xxtr;
                end
                disp('WTF');
            end
        end
    elseif and(xxtr(1)<-StdNoise,xxtr(end)>-StdNoise)
    % 1st Sample below & Final Sample above noise
        if or( numel(xxtr(xxtr>StdNoise))>numel(xxtr(xxtr<0)),...
            or(xxtr(1)>StdNoise,xxtr(end)>-StdNoise) )
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            disp('Ca2+ Possible Transient Rescued')
        else
            if numel(xlinc)>3
                PeakAmp=findpeaks(xxtr-xlinc);
                PeakAmp=PeakAmp(PeakAmp>StdNoise);
            else
                PeakAmp=[];
            end
            if isempty(PeakAmp) && sum(xlinc>xxtr)>numel(xxtr)/2
                % If line is above data, use smooth:
                xlinc=xxtr;
            end
        end
    else
        % if xxtr(end)>StdNoise
        %    mslope=0;
        %    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
        %    disp('Ca2+ Transient')
        % else
            % Otherwise Line
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            if numel(xlinc)>3
                PeakAmp=findpeaks(xxtr-xlinc);
                PeakAmp=PeakAmp(PeakAmp>StdNoise);
            else
                PeakAmp=[];
            end
            if isempty(PeakAmp) && sum(xlinc>xxtr)>numel(xxtr)/2
                % If line is above data, use smooth:
                xlinc=xxtr;
            end
        % end
    end
end