function sigma2(op::Array{Complex128,2})
    real(sum(diag(op^2) .- diag(op).^2))
end

function sigma2mult(ops::Array{Array{Complex128,2},1})
	sig2 = 0.0
	for o in ops
		sig2 += sigma2(o)
	end
	sig2
end

function sweepAJD!(ops::Array{Array{Complex128,2},1},U::Array{Complex128,2})
	NN = size(ops[1])[1] #Assume they are all of same size and square
	n = length(ops)
	for j = 1:NN
		for i = 1:(j-1)
			h = Array{Float64}(3)
			G = zeros(Float64,3,3)
			for opNum = 1:n
				h[1] = real(ops[opNum][i,i] - ops[opNum][j,j])
				h[2] = 2*real(ops[opNum][i,j])
				h[3] = 2*imag(ops[opNum][i,j])
				G += h*h'
			end
			F = eigfact(G)
			xx, yy, zz = F[:vectors][1,3] >= 0 ? F[:vectors][:,3] : -F[:vectors][:,3]
			rr::Float64 = sqrt(xx^2 + yy^2 + zz^2)
			c::Float64 = sqrt( (xx+rr)/(2*rr) )
			s::Complex128 = (yy - 1im*zz)/sqrt(2*rr*(xx+rr))
			R = Array{Complex128}(2,2)
			R[1,1] = c
			R[1,2] = conj(s)
			R[2,1] = -s
			R[2,2] = conj(c)
			for opNum = 1:n
				ops[opNum][:,[i,j]] = ops[opNum][:,[i,j]]*R'
				ops[opNum][[i,j],:] = R*ops[opNum][[i,j],:]
			end
			U[:,[i,j]] = U[:,[i,j]]*R'
		end
	end
end

function AJD(ops::Array{Array{Complex128,2},1};max_iter::Int64 = 500, tol::Float64 =1e-14)
	for A=ops
		NN, MM = size(ops])
		if NN!=MM
			error("At least one matrix is not square.")
		end
	end
	U = eye(Complex128,NN)
	newOps = copy(ops)
	postSweep::Float64 = sigma2mult(newOps)
	for ii = 1:max_iter
		preSweep::Float64 = postSweep
		sweepAJD!(newOps,U)
		postSweep = sigma2mult(newOps)
		if (preSweep - postSweep)<tol
			break
		elseif ii==max_iter
			info("Completed ",max_iter, " iterations without convergence: ",
				(preSweep - postSweep), " >= ", tol)
		end
	end
	newOps, U, postSweep
end

Upp(G::Array{Complex128,2}) = expm((G-G')/2)

function lowerSigma(eps::Float64,ops::Array{Array{Complex128,2},1},G::Array{Complex128,2})
    newops::Array{Array{Complex128,2},1} = copy(ops)
    c=1
    Utemp::Array{Complex128,2} = expm(eps*(G-G')/2)
    for XX in newops
        newops[c] = Utemp'*newops[c]*Utemp
        c+=1
    end
    sigma2mult(newops) - sigma2mult(ops)
end

function step_SteepestDescent!(ops::Array{Array{Complex128,2},1},U::Array{Complex128,2})
	NN = size(ops[1])[1] #Assume they are all of same size and square
	G = zeros(Complex128,NN,NN)
	n = length(ops)
	for opNum = 1:n
		G += 2*(ops[opNum]*spdiagm(diag(ops[opNum])) - spdiagm(diag(ops[opNum]))*ops[opNum])
	end

	optimalStep = optimize(eps -> lowerSigma(eps,ops,G),0,1)
	epsMin = 0.0
	if Optim.minimum(optimalStep) < 0.0
		epsMin = Optim.minimizer(optimalStep)
	end

	dU = expm(epsMin*(G-G')/2)

	for opNum = 1:n
		ops[opNum] = dU'*ops[opNum]*dU
	end

	U[:,:] = U*dU

	epsMin
end
