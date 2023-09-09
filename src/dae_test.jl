using LinearAlgebra

include("cross.jl")

N_b = 8;
N_DOF = N_b*6;



masses = rand(N_DOF)*50;

M = diagm(masses);

K = rand(N_DOF,N_DOF)*100; K = K*K';
# reduce the order of ODE
# M*ddu + C*du + K*u = 0
#
# u = u1 and du = u2;
# then 
#     du1 = 0 * u1 + 1 * u2;
#     du2 = M^(-1) * K * u1 + M^(-1) * C * u2

C = 0.1 *M +0.01*K;
A = [zeros(N_DOF, N_DOF) diagm(ones(N_DOF));
     inv(M)*K inv(M)*C];

function func(u)
    return A*u;
end

u0 = rand(2*N_DOF)*10;




t_step = 25/1e6;
using BenchmarkTools
@benchmark solve_cross(diagm(ones(2*N_DOF)), func, u0, t_step, (0, 0.1))
@benchmark solve_cross_linear(diagm(ones(2*N_DOF)), A, u0, t_step, (0, 0.1))
# t, u = solve_cross(diagm(ones(2*N_DOF)), func, u0, t_step, (0, 1));

plot(t, u[1, :])

# plot3D(u[1,:], u[2,:], u[3,:])


