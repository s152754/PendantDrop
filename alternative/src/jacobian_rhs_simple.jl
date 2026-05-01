function jacobian_rhs_simple(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Matrix,Vector}
    
    D = vars_num.D0;
    w = vars_num.w0;
    N = vars_num.N;

    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;

    g = params_phys.grav;
    sigma = params_phys.sigma;
    rho = params_phys.deltarho;
    V = params_phys.volume0;
    
    # initialize some variables 
    Z = zeros(N,N);            # matrix filled with zeros
    IDL = [1 zeros(1,N-1)]; # line with single one and rest zeros
    ZL = zeros(1,N);         # line completely filled with zeros 
    b = ones(3*N+2,1); # solution vector and right hand side
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    # determine r from psi
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    
    # determine z from psi 
    A22 = C*D;
    A23 = Diagonal(-cos.(psi));
    A24 = D*z;
    b2 = -(C*D*z-sin.(psi));
    
    # determine psi from Laplace law
    A31 = -sigma*Diagonal(sin.(psi)./r.^2);
    A32 = g*rho*Diagonal(ones(N));
    A33 = C*sigma*D + sigma*Diagonal(cos.(psi)./r);
    A34 = sigma*(D*psi);
    A35 = -ones(N,1);
    b3 = p0.-g*rho*z-sigma*(C*D*psi+sin.(psi)./r);
    
    # impose the needle radius as a BC (imposes the domain length)
    A41 = reverse(IDL);
    b4 = (params_phys.rneedle-r[end]);
    
    # determine pressure - use volume
    A51 = pi*2*(w.*r.*sin.(psi))';
    A53 = pi*(w.*r.^2 .*cos.(psi))';
    A54 = -V;
    b5 = -(pi*w'*(r.^2 .*sin.(psi)).-C*V);
    
    # boundary condition r(0) = 0
    A11[1,:] = IDL;
    A13[1,:] = ZL; 
    A14[1] = 0;
    b1[1] = -r[1];
    
    # boundary condition z(s0) = 0
    A22[1,:] = reverse(IDL); 
    A23[1,:] = ZL; 
    A24[1] = 0;
    b2[1] = -z[end];
    
    # boundary condition phi(0) = 0
    A31[1,:] = ZL;
    A32[1,:] = ZL;

    A33[1,:] = IDL; 
 
    A34[1,:] .= 0; 
    A35[1,:] .= 0;
    b3[1] = -psi[1];

    
    # assemble matrices
    Z1 = zeros(N,1);
     
    A = [[A11   Z A13 A14  Z1];
       [  Z A22 A23 A24  Z1];
       [A31 A32 A33 A34 A35];
       [A41  ZL  ZL   0   0];
       [A51 Z1' A53 A54   0]];
     
    b = [b1;b2;b3;b4;b5];

    return (A,b)
end