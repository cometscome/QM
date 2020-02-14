N=1000
a = 0.01
dx = N/6
ξ = dx*a
center = (N+1)/2
x0 = center*a

function calc_vq(q,ξ,V0)
    vq = sqrt(π*ξ^2)*exp(-q^2*ξ^2/4)
    return vq
end

function calc_Vkkp(k,kp,ξ,x0,V0)
    q1 = k - kp
    vq1 = calc_vq(q1,ξ,V0)
    q2 = k + kp
    vq2 = calc_vq(q2,ξ,V0)
    Vkkp = 2*V0*(cos(q1*x0)*vq1 - cos(q2*x0)*vq2)
    return Vkkp
end

function make_Hk(N,a,V0)
    mat_Hk = zeros(Float64,N,N)
    dx = N/6
    ξ = dx*a
    center = (N+1)/2
    x0 = center*a
    L = (N+1)*a
    for n in 1:N
        k = n*π/L
        for np in 1:N
            v = 0.0            
            if n == np
                v = k^2
            end
            kp = np*π/L
            Vkkp = calc_Vkkp(k,kp,ξ,x0,V0) 
            v += Vkkp*(1/2L)
            mat_Hk[n,np]= v
        end
    end
    return mat_Hk
end

using LinearAlgebra #対角化のルーチンeigenを呼ぶ準備 v0.7以降必要

V0 = 0.0
mat_H = make_Hk(N,a,V0)
ε,ψ = eigen(mat_H)
println("Potential = ",V0," Minimum eigenvalue = ",ε[1])

a = 0.01
L = (N+1)*a
ε1 = π^2/L^2
println(ε1)

using Plots
gr()

V0 = 1.0
N = 1000
a = 0.01
mat_H = make_Hk(N,a,V0)
ep,psi = eigen(mat_H)
println(ep[1])
a = 0.01
rp = zeros(Float64,N)
L = a*N
for i in 1:N
    xi = a*i
    for ik in  1:N
        k = π*ik/L
        rp[i] += psi[ik,1]*sin(k*xi)
    end
end


C = sum(dot(rp[1:N],rp[1:N]))
rp = rp/sqrt(C)

xin = []
for i in 1:N
    push!(xin,i)
end
plot(xin,rp)
savefig("fig03_1.png")



function make_H1dv(N,a,V0)
    mat_H = zeros(Float64,N,N)
    vec_V = calc_V(N,V0)
        
    for i in 1:N
        for dx in -1:1
            j = i + dx
            v = 0.0
            if dx == 0
                v = (2/a^2 + vec_V[i])
            elseif dx == 1
                v = -1/a^2
            elseif dx == -1
                v = -1/a^2
            end
            
            if 1 <= j <= N
                mat_H[i,j] = v
            end
            
        end
        
    end
    
    
    return mat_H
end

function calc_V(N,V0)
    vec_V = zeros(Float64,N)
    dx = N/6
    center = (N+1)/2
    for i in 1:N
        vec_V[i] = V0*exp(-(i-center)^2/(dx^2))
    end
    return vec_V
end

N = 1000
a = 0.01
V0=1.0
mat_Hx = make_H1dv(N,a,V0)
εx,ψx = eigen(mat_Hx)
plot(xin,[rp,ψx[1:N,1]],label=["momentum based method" "x-based method"])

savefig("fig03_2.png")

a = 0.01
N=1000

function gs1(a,N)
    minimums =[]
    for v in 1:10
        V0 = v*0.5
        mat_H = make_Hk(N,a,V0)
        ε,ψ = eigen(mat_H)
        push!(minimums,ε[1])
        println("Potential = ",V0," Minimum eigenvalue = ",ε[1])
    end  
    return minimums
end
minimums = gs1(a,N)

a = 0.01
N=1000
function gs2()
    minimums_x =[]
    for v in 1:10
        V0 = v*0.5
        mat_Hx = make_H1dv(N,a,V0)
        εx,ψx = eigen(mat_Hx)
        push!(minimums_x,εx[1])
        println("Potential = ",V0," Minimum eigenvalue = ",εx[1])
    end 
    return  minimums_x
end
minimums_x = gs2()

function gs3()
    potentials = []
    for v in 1:10
        V0 = v*0.5
        push!(potentials,V0)
    end
    return potentials
end
potentials = gs3()

plot(potentials,[minimums,minimums_x],label=["k-based method" "x-based method"])
savefig("fig03_3.png")


V0 = 40.0
N = 1000
a = 0.01
mat_H = make_Hk(N,a,V0)
ep,psi = eigen(mat_H)
println("momentum-based method: ",ep[1])
a = 0.01
rp = zeros(Float64,N)
L = a*N
for i in 1:N
    xi = a*i
    for ik in  1:N
        k = π*ik/L
        rp[i] += psi[ik,1]*sin(k*xi)
    end
end


C = sum(dot(rp[1:N],rp[1:N]))
rp = rp/sqrt(C)

xin = []
for i in 1:N
    push!(xin,i)
end

mat_Hx = make_H1dv(N,a,V0)
εx,ψx = eigen(mat_Hx)
println("x-based method ",εx[1])
plot(xin,[rp,ψx[1:N,1]],label=["momentum based method" "x-based method"])
savefig("fig03_4.png")

V0 = 40.0
N = 2000
a = 0.005
mat_H = make_Hk(N,a,V0)
ep,psi = eigen(mat_H)
println("momentum-based method: ",ep[1])
a = 0.01
rp = zeros(Float64,N)
L = a*N
for i in 1:N
    xi = a*i
    for ik in  1:N
        k = π*ik/L
        rp[i] += psi[ik,1]*sin(k*xi)
    end
end


C = sum(dot(rp[1:N],rp[1:N]))
rp = rp/sqrt(C)

xin = []
for i in 1:N
    push!(xin,i)
end

mat_Hx = make_H1dv(N,a,V0)
εx,ψx = eigen(mat_Hx)
println("x-based method ",εx[1])
plot(xin,[rp,ψx[1:N,1]],label=["momentum based method" "x-based method"])
savefig("fig03_5.png")


a = 0.01
N=1000
function gs4()
    minimums =[]
    V0 = 40.0
    for nn in 1:10
        N = 250*nn
        a = 10.0/N
        mat_H = make_Hk(N,a,V0)
        ε,ψ = eigen(mat_H)
        push!(minimums,ε[1])
        println("Number = ",N," Minimum eigenvalue = ",ε[1])
    end 
    return minimums,ψ
end
minimums,ψ = gs4()

a = 0.01
N=1000
function gs5()
    minimums_x =[]
    for nn in 1:10
        N = 250*nn
        a = 10.0/N
        mat_Hx = make_H1dv(N,a,V0)
        εx,ψx = eigen(mat_Hx)
        push!(minimums_x,εx[1])
        println("Number = ",N," Minimum eigenvalue = ",εx[1])
    end  
    return minimums_x,ψx
end
minimums_x,ψx = gs5()

function gs6()
    numbers = []
    for nn in 1:10
        N = 250*nn
        push!(numbers,N)
    end
    return numbers
end
numbers = gs6()

plot(numbers,[minimums,minimums_x],label=["k-based method" "x-based method"])
savefig("fig03_6.png")


N = 2500
V0 = 40.0
a = 10.0/N
mat_H = make_Hk(N,a,V0)
ε,ψ = eigen(mat_H)

mat_Hx = make_H1dv(N,a,V0)
εx,ψx = eigen(mat_Hx)

rp = zeros(Float64,N)
L = a*N
for i in 1:N
    xi = a*i
    for ik in  1:N
        k = π*ik/L
        rp[i] += ψ[ik,1]*sin(k*xi)
    end
end


C = sum(dot(rp[1:N],rp[1:N]))
rp = rp/sqrt(C)

xin = []
for i in 1:N
    push!(xin,i)
end

println("x-based method ",εx[1])
plot(xin,[-rp,ψx[1:N,1]],label=["momentum based method" "x-based method"])
savefig("fig03_7.png")