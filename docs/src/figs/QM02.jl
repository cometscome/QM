function make_H1d(N,a)
    mat_H = zeros(Float64,N,N)
    vec_V = zeros(Float64,N)
        
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

using LinearAlgebra #対角化のルーチンeigenを呼ぶ準備 v0.7以降必要
N = 1000
a = 0.01
mat_H = make_H1d(N,a)
ε,ψ = eigen(mat_H)

println(ε[1])

is = 1
println(ψ[:,is])

using Plots
gr()

xaxis = Int64[]
for i in 1:N
    push!(xaxis,i)
end

is = 1
plot(xaxis[1:N],ψ[1:N,is],label="Eigenfunction") 
savefig("fig02_1.png")


aψ = zeros(Float64,N)
for i in 1:N
    xi = i*a
    aψ[i] = sin(xi*π/((N+1)*a))
end

is = 1
plot(xaxis[1:N],[ψ[1:N,is],aψ[1:N]],label=["Numerical result" "Analytical result"]) 
savefig("fig02_2.png")


aψ = zeros(Float64,N)
for i in 1:N
    xi = i*a
    aψ[i] = sin(xi*π/((N+1)*a))
end
C = sum(dot(aψ[1:N],aψ[1:N]))
aψ = aψ/sqrt(C)
is = 1
plot(xaxis[1:N],[ψ[1:N,is],aψ[1:N]],label=["Numerical result" "Analytical result"]) 
savefig("fig02_3.png")

plot(xaxis[1:N],ψ[1:N,1:6],label=["1st" "2nd" "3rd" "4th" "5th" "6th"]) 
savefig("fig02_4.png")

xaxis = Int64[]
for i in 1:N
    push!(xaxis,i)
end

function calc_V(N,V0)
    vec_V = zeros(Float64,N)
    dx = N/6
    for i in 1:N
        if N/2 - dx <= i <= N/2 + dx
            vec_V[i] = V0
        end
    end
    return vec_V
end

V0=1.0
vec_V = calc_V(N,V0)
plot(xaxis[1:N],vec_V[1:N],label="Potential") 
savefig("fig02_5.png")

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

N = 1000
a = 0.01
V0 = 1.0
mat_H = make_H1dv(N,a,V0)
ε,ψ = eigen(mat_H)

println(ε[1])
is = 1
plot(xaxis[1:N],ψ[1:N,is],label="Eigenfunction") 


savefig("fig02_6.png")


N = 1000
a = 0.01
function gs1(N,a)
    groundstates = []
    labels = []
    for v in 1:10
        V0 = v*0.5
        mat_H = make_H1dv(N,a,V0)
        ε,ψ = eigen(mat_H)
        println("Potential = ",V0," Minimum eigenvalue = ",ε[1])
        push!(groundstates,ψ[:,1])
        if v == 1
            labels = [string(V0)]
        else
            labels = [labels string(V0)]
        end
        #push!(labels,string(V0))
    end   
    return groundstates,labels
end
groundstates,labels = gs1(N,a)

plot(xaxis[1:N],groundstates,label=labels) 
savefig("fig02_7.png")

function calc_V(N,V0)
    vec_V = zeros(Float64,N)
    dx = N/6
    center = (N+1)/2
    for i in 1:N
        if center - dx <= i <= center + dx
            vec_V[i] = V0
        end
    end
    return vec_V
end


N = 1000
a = 0.01
V0 = 5.0
mat_H = make_H1dv(N,a,V0)
ε,ψ = eigen(mat_H)

println(ε[1])
is = 1
plot(xaxis[1:N],ψ[1:N,is],label="Eigenfunction") 
savefig("fig02_8.png")



function calc_V(N,V0)
    vec_V = zeros(Float64,N)
    dx = N/6
    center = (N)/2
    for i in 1:N
        if center - dx <= i <= center + dx
            vec_V[i] = V0
        end
    end
    return vec_V
end


N = 1000
a = 0.01
V0 = 5.0
mat_H = make_H1dv(N,a,V0)
ε,ψ = eigen(mat_H)
println(ε[1])
plot(xaxis[1:N],ψ[1:N,1],label="a=0.01") 
savefig("fig02_9.png")

N = 1000*2
xaxis = Int64[]
for i in 1:N
    push!(xaxis,i)
end


a = 0.01/2
V0 = 5.0
mat_H = make_H1dv(N,a,V0)
println(ε[1])
ε,ψ = eigen(mat_H)
plot(xaxis[1:N],ψ[1:N,1],label="a=0.005") 
savefig("fig02_10.png")


function calc_V(N,V0)
    vec_V = zeros(Float64,N)
    dx = N/6
    center = (N)/2
    for i in 1:N
        vec_V[i] = V0*exp(-(i-center)^2/(dx^2))
    end
    return vec_V
end

N = 1000
xaxis = Int64[]
for i in 1:N
    push!(xaxis,i)
end

V0=1.0
vec_V = calc_V(N,V0)
plot(xaxis[1:N],vec_V[1:N],label="Potential") 
savefig("fig02_11.png")

N = 1000
a = 0.01
function gs2(N,a)
    groundstates = []
    labels = []
    for v in 1:10
        V0 = v*0.5
        mat_H = make_H1dv(N,a,V0)
        ε,ψ = eigen(mat_H)
        println("Potential = ",V0," Minimum eigenvalue = ",ε[1])
        push!(groundstates,ψ[:,1])
        if v == 1
            labels = [string(V0)]
        else
            labels = [labels string(V0)]
        end
    end    
    return groundstates,labels
end
groundstates,labels = gs2(N,a)

plot(xaxis[1:N],groundstates,label=labels) 
savefig("fig02_12.png")

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
function gs3(N,a)
    groundstates = []
    labels = []
    for v in 1:10
        V0 = v*0.5
        mat_H = make_H1dv(N,a,V0)
        ε,ψ = eigen(mat_H)
        println("Potential = ",V0," Minimum eigenvalue = ",ε[1])
        push!(groundstates,ψ[:,1])
        if v == 1
            labels = [string(V0)]
        else
            labels = [labels string(V0)]
        end
    end    
    return groundstates,labels
end
groundstates,labels = gs3(N,a)

plot(xaxis[1:N],groundstates,label=labels) 
savefig("fig02_13.png")