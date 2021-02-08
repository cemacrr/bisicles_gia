def post(c_one,mu_coef,u,v,uo,vo,uc,*etc):

    um = (u*u + v*v)**0.5
    uo = (uo*uo + vo*vo)**0.5
    c_third = c_one * (1+um)**(2.0/3.0)

    return c_one, c_third, mu_coef, um, uo, uc
