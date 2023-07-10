#------------------------------------------------------------------------------------------------------------------------------
# Module Harmonic
#------------------------------------------------------------------------------------------------------------------------------
# Package Parameters
# julia> using Pkg
# julia> Pkg.add("Parameters")
#
# Package Plots
# julia> Pkg.add("Plots")
#
# Packages Polynomials and Special Polynomials
# julia> import Pkg
# julia> Pkg.add("Polynomials")
# julia> Pkg.add("SpecialPolynomials")

module Harmonic

using LinearAlgebra
using Plots
using Polynomials, SpecialPolynomials

using ..constants

# How to find the path of a package in Julia?
# https://stackoverflow.com/questions/61440940/how-to-find-the-path-of-a-package-in-julia

# Function module file looks up the location of a module from the method table
module_file(modu) = String(first(methods(getfield(modu, :eval))).file)

println("File Harmonic.jl")
println("Using constants from file: ", module_file(constants))

export mu_co2, omegacm1, omega, r_co
export V_harm, T_harm, Hop_grid, H_op, vib_freq_CO2, convergence_test, WF_harm_osc 

# Define constants
const mu_co2 = 2*constants.mO::Float64
const omegacm1 = 1333.0::Float64
const omega = omegacm1*constants.cm1::Float64
const r_co = 2.1958618::Float64 # equilibrium bond length in bohrs

# Define function V_harm(Q,omg_nu, mu_q)
# We use the experimental frequency to construct the harmonic potential. The experimental frequency is taken from NIST. 
# Most of websites show slightly different value for the frequency. It might be because of the chosen basis and parameters 
# in the quantum chemistry packages they are using. However, we use NIST value as the main reference.
function V_harm(Q, omg_nu::Float64, mu_q::Float64)
    fac::Float64 = 0.5*mu_q*(omg_nu^2)
    return fac*(Q).^2
end

# Define function T_harm(QQ,mu_q)
# Compute the kinetic energy matrix over the grid of QQ
# Define the kinetic energy matrix
# Compute the kinetic energy matrix over the grid of QQ
# Define the kinetic energy matrix
function T_harm(QQ::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}, mu_q::Float64)
    
    N = length(QQ)
    h = QQ[2]-QQ[1]
    
    T0 = 2 .* ones(1,N)
    T1 = -1 .* ones(1,N-1)
    Tm1 = -1 .* ones(1,N-1)
    
    Tu = Tridiagonal(vec(Tm1),vec(T0),vec(T1))
        
    ## Apply the boundary conditions
    Tu[1,1] = 0
    Tu[1,2] = 0
    Tu[2,1] = 0
    Tu[N-1,N-2] = 0
    Tu[N-2,N-1] = 0
    Tu[N-1,N-1] = 0
    
    T = Tu/(2*mu_q*(h^2))
    
    return T
end


# Define function Hop_grid(QQ,omega,mu)
function Hop_grid(QQ::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},omega::Float64,mu::Float64)   

    # Initiate the matrix
    QQn = length(QQ)

    # The potential energy surface over the grid of QQ
    VV::Matrix{Float64} = diagm([V_harm(i,omega,mu) for i in QQ])
    V::Diagonal{Float64, Vector{Float64}} = Diagonal(VV)

    # Kinetic Energy matrix
    T::Tridiagonal{Float64, Vector{Float64}} = T_harm(QQ,mu)
    
    # Construct the Hamiltonian 
    H = T+V
    
    return H
end

# Define function vibrational_freq_CO2() to compute the vibrational frequency of  𝐶𝑂2
# Compute the vibrational frequency of the CO2 molecule on the symmetric stretch normal mode
function vib_freq_CO2(H_op::Tridiagonal{Float64, Vector{Float64}}, omega::Float64, len::Int64)
    
    Evalues, Evecs = eigen(H_op)
    
    # make sure that the eigenvalues are sorted in ascending order
    Evalues = sort(Evalues)
    
    # convert the eigenvalues from Eh to cm^-1
    Ecm1 = Evalues./constants.cm1
    
    # compute the error
    error = ((Ecm1[2]-Ecm1[1])-omegacm1)/omegacm1
    
    return Evalues, Ecm1, error
end

# Define function convergence test(Npoints, omega, mu) to test the convergence based on the number of points in the grid of QQ
function convergence_test(Npoints::Vector{Int64}, omega::Float64, mu::Float64)

    len = size(Npoints)[1]

    # Define an array for the errors
    error = zeros(len)

    i = 1
    for l in 1:len
        Qvec = range(start=first(Npoints), stop=last(Npoints), step=1/l)
    
        # Call the Hamiltonian function
        H = Hop_grid(Qvec, omega, mu)
    
        # Calculate the vibrational frequency and the error
        Evalues, Ecm1, error[i] = vib_freq_CO2(H, omega, length(Qvec))
        
        i += 1
    end

    return error
end

## Define harmonic_oscillator_wavefunction(x, n, m, omega):
# Hermite polynomials, https://en.wikipedia.org/wiki/Hermite_polynomials
# Julia Special Polynomials, https://docs.juliahub.com/SpecialPolynomials/LrhA0/0.1.0/

function WF_harm_osc(x::Float64, n::Int, mu::Float64, omega::Float64)
    # Inputs:
    # x: position
    # n: quantum number
    # mu: reduced mass
    # omega: vibrational frequency
    # Returns:
    # psi: the harmonic wavefunction
    prefactor = sqrt(1.0/(2^n * factorial(n))) * ((mu*omega) / pi)^0.25
    exponent = -0.5 * mu * omega * (x)^2
    herm = basis(Hermite, n)
    return prefactor * herm(sqrt(mu * omega) * x) * exp(exponent)
end

end # module Harmonic