function D = Diffusion(Re, F, parameters)

    a = parameters{44};
    b = parameters{45};
    c = parameters{46};

    D =  a -  b * Re + c  * F * 10^5;
    
end