function [P0, P1, Pm1, Pt] = probabilitycalculatorv1(H, initial, t)

    digits(7)
    %% eigenvectors/values of hamiltonian H
    [vvecs, vvals] = eig(H);
    %% define unit vectors
    e1 = [1;0;0];
    e2 = [0;1;0];
    e3 = [0;0;1];
    %% solve to find e2 in the energy basis
    E123 = linsolve(vvecs, initial);
    %% create components of psi(t)
    E1sol = E123(1)*exp(-1i*vvals(1,1)*t*(1/0.6582));
    E2sol = E123(2)*exp(-1i*vvals(2,2)*t*(1/0.6582));
    E3sol = E123(3)*exp(-1i*vvals(3,3)*t*(1/0.6582));
    %% decompose each into its spin-basis components
    S1sol = [E1sol*vvecs(1,1), E1sol*vvecs(2,1), E1sol*vvecs(3,1)];
    S2sol = [E2sol*vvecs(1,2), E2sol*vvecs(2,2), E2sol*vvecs(3,2)];
    S3sol = [E3sol*vvecs(1,3), E3sol*vvecs(2,3), E3sol*vvecs(3,3)];
    %% add together to obtain each spin component
    Ssol = S1sol + S2sol + S3sol;
    %% find probabilities of each state, P0 is probability it's still in the ground state
    P0 = dot(e2,Ssol)*conj(dot(e2, Ssol))
    P1 = (abs(dot(e1,Ssol)))^2
    Pm1 = (abs(dot(e3,Ssol)))^2
    %% make sure that all probabilities sum to 1!
    Pt = P0 + P1 + Pm1

end
