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

using Plots
gr()

using LinearAlgebra #対角化のルーチンeigenを呼ぶ準備 v0.7以降必要

N = 1000
a = 0.01
mat_H = make_H1d(N,a)
ε,ψ = eigen(mat_H)

integers = Int64[]
for i in 1:N
    push!(integers,i)
end

plot(integers[1:N],ε[1:N],label="1D") 
savefig("fig01.png")

sn = 200
plot(integers[1:sn],ε[1:sn],label="1D",marker=:circle) 
savefig("fig02.png")


εa = Float64[]
for n in 1:N
    push!(εa,n^2*π^2/(a*(N+1))^2)
end

plot(integers[1:sn],[ε[1:sn],εa[1:sn]],label=["Numerical result" "Analytical result"],marker=:circle) 
savefig("fig03.png")

sn = 1000
plot(integers[1:sn],[ε[1:sn],εa[1:sn]],label=["Numerical result" "Analytical result"],marker=:circle)
savefig("fig04.png")

N = 1000
a = 0.1
mat_H = make_H1d(N,a)
ε,ψ = eigen(mat_H)

integers = Int64[]
for i in 1:N
    push!(integers,i)
end

sn = 400
εa = Float64[]
for n in 1:N
    push!(εa,n^2*π^2/(a*(N+1))^2)
end

plot(integers[1:sn],[ε[1:sn],εa[1:sn]],label=["Numerical result" "Analytical result"],marker=:circle) 

savefig("fig05.png")