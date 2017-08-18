function [C,angles] = prinAngles(A,B, method)

if strcmp(method,'fast')
    % quick algorithm for calculating the principal angles
    % between the two subspaces spanned by the columns of
    % A and B when angles are larger than 10^(-6).

    [Qa,Ra] = qr(A,0);
    [Qb,Rb] = qr(B,0);
    C = svd(Qa'*Qb,0);
    angles = acos(C);

else

    % a modified algorithm for calculating the principal angles %
    % between the two subspaces spanned by the columns of       %
    % A and B. Good for small (<10^(-6)) and large angles.   %
    
    [Qa,Ra] = qr(A,0);
    [Qb,Rb] = qr(B,0);
    C = svd((Qa')*Qb,0);
    rkA = rank(Qa);
    rkB = rank(Qb);
    if rkA >= rkB
    B = Qb - Qa*(Qa'*Qb);
    else
    B = Qa - Qb*(Qb'*Qa);
    end
    S = svd(B,0);
    S = sort(S);
    for i = 1:min(rkA,rkB)
        if (C(i))^2 < 0.5
            angles(i) = acos(C(i));
        elseif (S(i))^2 <= 0.5
            angles(i) = asin(S(i));
        end
    end
    angles=angles';

end