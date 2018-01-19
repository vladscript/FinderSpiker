%% Function to make Row Vector any vector
% Input: vector A
% Output: vector B
function B=makerowvector(A)
if isvector(A)
    if isrow(A)
        B=A;
    else
        B=A';
    end
else
    disp(' > Just Tranpose Matrix < ')
    B=A';
end