function xn = myCramer_PGBC(A, b)
% input b should be a column vector

lenA = length(A);
xn = zeros(lenA,1);
detA = det(A);
for i = 1:lenA
    Ai = [A(:,1:i-1) b A(:,i+1:lenA)];
    detAi = det(Ai);
    xn(i) = detAi/detA;
end

% Presented by Pio Gabrielle B Calderon