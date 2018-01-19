%% Wavelet Anlysis for Calcium Imaging Signals
% Input
%   xd:         Signal
%   wname:      Wavelet Name
%   level:      Level of Analysis
% Output
%   xdenoised:  Denoised Signal
function sigden=waveletdenoise(xd,wname,level)
% level=1;
% wname='sym8';
[coefs,longs] = wavedec(xd,level,wname);
% siz = size(coefs);
% Error here:-------------\/

[thr,sorh,keepapp] = ddencmp('den','wv',xd);

for i=1:level
    thrParams{i}=[size(xd),thr];
end
% thrParams = utthrset_cmd(coefs,longs);
first = cumsum(longs)+1;
first = first(end-2:-1:1);
tmp   = longs(end-1:-1:2);
last  = first+tmp-1;
for k = 1:level
    thr_par = thrParams{k};
    if ~isempty(thr_par)
        cfs = coefs(first(k):last(k));
        nbCFS = longs(end-k);
        NB_int = size(thr_par,1);
        x = [thr_par(:,1) ; thr_par(NB_int,2)];
        alf = (nbCFS-1)/(x(end)-x(1));
        bet = 1 - alf*x(1);
        x = round(alf*x+bet);
        x(x<1) = 1;
        x(x>nbCFS) = nbCFS;
        thr = thr_par(:,3);
        for j = 1:NB_int
            if j==1 ,
                d_beg = 0;
            else
                d_beg = 1;
            end
            j_beg = x(j)+d_beg;
            j_end = x(j+1);
            j_ind = (j_beg:j_end);
            cfs(j_ind) = wthresh(cfs(j_ind),sorh,thr(j));
        end
        coefs(first(k):last(k)) = cfs;
    end
end
sigden = waverec(coefs,longs,wname);