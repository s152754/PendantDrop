using LinearAlgebra
using Interpolations


include("PD_functions.jl")
include("PD_structs.jl")
include("utils.jl")
include("plotting.jl")

function PD_fn(params; rneedle=1., deltarho=1.26e-6, grav=9.81e3, N=100, nMaxIter=100, epsilon=1e-12)
    # Unpack parameters
    sigma = params[1]  # Surface tension (mN/m)
    volume0 = params[2]  # Prescribed volume (mm^3)

    # Assign structs containing the parameters you would want to solve for
    params_phys = ParamsPhys(sigma=sigma, grav=grav, rneedle=rneedle, volume0=volume0, deltarho=deltarho)
    params_num = ParamsNum(N=N, nMaxIter=nMaxIter, epsilon=epsilon)

    # Solve for the droplet shape (Young-Laplace)
    Ngrid = 40
    shape_guess = guess_shape(params_phys, Ngrid)
    vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=false, shape_guess = shape_guess)
    volume, area = calculate_volume_area(vars_sol, vars_num; verbose=false);
    kappas, kappap = find_curvature(vars_sol, vars_num);

    return vars_sol.r, vars_sol.z, volume, area, kappas, kappap
end