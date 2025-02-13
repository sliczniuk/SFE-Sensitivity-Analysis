function [T] = reconstruct_T_polynomial_approximation(H, P)

    T = 190.1 + 0.04924.*H + 2.651.*P + 0.006511.*H.^2 + 0.01881.*H.*P + -0.001048.*P.^2 + 2.098e-05.*H.^3 + 4.062e-05.*H.^2.*P + 2.552e-06.*H.*P.^2 + 2.128e-06.*P.^3;

end