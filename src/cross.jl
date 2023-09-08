# With give equation 
# M * du/dt = F(u)
# u(0) = u_0
# the solution can be found in form 
# u_(i+1) = u_i + dt * Re(k)
# where k - is a solution of the linear equation (M - 0.5(j+1) * J_F * dt) * k = F(u_i + dt/2)
# dt is a time step, j = sqrt(-1) - imaginary number

using ForwardDiff
using LinearAlgebra
using PyPlot


function solve_cross(M, F, u0, dt, time_span)
    # construct the jacoby matrix of form
    
    J(x) = ForwardDiff.jacobian(F, x);
    t = time_span[1]:dt:time_span[2];
    u = Matrix{Float64}(undef, length(u0), length(t));
    u[:,1] = u0;

    q = 0.5*(1+im);
    
    for i in 2:length(t)
        A = M - q*dt*J(u[:,i-1])
        k = A\F(u[:,i-1]);
        u[:,i] = u[:,i-1] + dt*real(k);
    end
    return t,u
end

function lorenz(u)
    du = similar(u);
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return du
end

function run_example(D3)

    u0 = [1.0, 0.0, 0.0];

    M = [1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0];

    t, u = solve_cross(M, lorenz, u0, 0.001, (0, 100));
    if !(D3)

        plot(t, u[1, :])
    else
        plot3D(u[1,:], u[2,:], u[3,:])
    end
end
