function update_numerical_grid(
    vars_sol::VarsSol,
    vars_num::VarsNum)::VarsNum

    # the integration and differentation matrices in the solution state
    vars_num.ws = vars_num.w0/vars_sol.C;
    vars_num.Ds = vars_sol.C*vars_num.D0; 
    vars_num.s = vars_num.s0/vars_sol.C;
    #vars_num.wsmat = vars_num.wmat/vars_sol.C;
    vars_num.C = vars_sol.C;

    return vars_num
end