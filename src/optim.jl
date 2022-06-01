function ConstantPreprocess(A::SparseMatrixCSC{Int64,Int64},d::Array{Int64,2},s::Float64,fast::Int64)

    N::Int64 = length(d); 
    v = rand(Uniform(-1,1),N,1);
    v = v/norm(v,2);
    f::Float64 = 0.;
    n::Float64 = 0.;
    
    if fast==1
        for it = 1:30
    
            v = A*v .- (d/s)*(d'*v);
            n = norm(v,2);
            v = v/n;
        end
        f = n;
    else
        lam = eigs(A-(d/s)*d',nev=1);
        f = real(lam[1][1]);
    end
    
    return f
    end
    


# Input: 
#
# A : sparse adjacency matrix N x N
# H0 : feasible initial matrix N x r0
# n_it : maximum number of iterations
# t : tolerance (used for assessing the stationarity of the sequence of iterates).
# fast: if ==1 this option approximates the largest eigenvalue of the modularity matrix (for large matrices)
# fast: if ==0 this option calculates the largest eigenvalue of the modularity matrix (for small matrices)
#
# Output:
#
# X : a feasible matrix N x r0
#########################################################

function AccProjPowMethodPreProcessed(A::SparseMatrixCSC{Int64,Int64},X::Array{Float64,2},n_it::Int64,t::Float64,fast::Int64)
   
d::Array{Int64,2} = sum(A,dims=2);
s::Float64 = sum(d);

# Initialization
N::Int64 = size(X,1);

f = vec(1 .+ 2*d -(d.^2)/s);

i::Int64 = 1;
diff::Float64 = 10.;

Y0::Array{Float64,2} = zeros(size(X));
Y::Array{Float64,2} = zeros(size(X));

o0::Float64 = 1.;
o::Float64 = 1.;
r::Float64 = 0.;

while (i<=n_it)  && (diff >t || i<4)	

   # Matrix product
   Y = A*X-(d/s)*(d'*X)+Diagonal(f)*X;
   # In order to avoid calculating twice the same products
   # Since X is feasible, we calculate the objective value AFTER the first matrix product.
   o=(tr(X'*Y))/s;  

   #Difference between consecutive objectives
   diff = abs.(o-o0)/o0; 

   #Update with momentum
   r = (i-1.)/(i+2.);
   X = (1 +r)*Y;
   X = X - r*Y0;

   #Project on the sphere
   X = ProjectSphere(X); #X = X./sqrt.(sum(abs2, X, 2));

   #Update objective and momentum
   o0 = o;
   Y0 = Y;

#	if mod(i,50)==0
	#	print("iteration: ")
	#	print(i)
	#	print("\n")
#	end

   i+=1;


end



# Output
if diff<t
  	@printf "The iteration has become stationary after %d iterations\n" i;
else
	@printf "The iteration did not converge after %d iterations\n" n_it;
	@printf "The relative difference between the last objective values %f \n" diff;
end


return X;
end
################################################################
