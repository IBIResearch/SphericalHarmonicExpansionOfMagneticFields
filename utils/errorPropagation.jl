"""
    coeffs, coeffsError, expansion, func = errorPropagation(coeffs::Array{SphericalHarmonicCoefficients}, 
                          coords::Array{T,2}, R::T, center::Vector{T}, fieldsError::Union{Array{T,2},Array{T,3}}, 
                          x::Variable, y::Variable, z::Variable; 
                          transGauss::Bool=true, transFFP::Bool=true, ffp=[0.0,0.0,0.0]) where T <: Real
"""
function errorPropagation(coeffs::Array{SphericalHarmonicCoefficients}, 
			  coords::Array{T,2}, R::T, center::Vector{T}, fieldsError::Union{Array{T,2},Array{T,3}}, 
			  x::Variable, y::Variable, z::Variable; 
			  transGauss::Bool=true, transFFP::Bool=true, ffp=[0.0,0.0,0.0]) where T <: Real

    # Translation in the Hall sensor
    v = zeros(3,3)
    v[:,1] = [0.00208, 0.0, -0.0018]
    v[:,2] = [0.0, 0.00208, -0.0018]
    v[:,3] = [0.0, 0.0, -0.0018]

    # Coordinate transformation Hall sensor -> Bruker Scanner:
    M = [0 0 1; 0 1 0; -1 0 0]
    v = M*v

    coeffsError = deepcopy(coeffs)
    func= Array{Function}(undef,3,size(coeffs,2))
    expansion = Array{Polynomial}(undef,3,size(coeffs,2))

    # rescale coordinates to t-design on unit sphere
    coords = coords .- center
    coords *= 1/R
    for j = 1:3
        for c = 1:size(coeffs,2)
            # Error of each coefficient
            coeffsError[j,c] = SphericalHarmonicExpansions.errorSphericalQuadrature(fieldsError[j,:,c],coords',coeffs[1,1].L);
            coeffsError[j,c].R = R

            # Scale error
            normalize!(coeffsError[j,c],R)
            if coeffs[j,c].solid
                coeffsError[j,c].solid = false
                solid!(coeffsError[j,c])
            end

            # Translation of Hall sensors
            if transGauss
                coeffsError[j,c] = SphericalHarmonicExpansions.errorTranslation(coeffsError[j,c],v[:,j])
            end

            # Translation to the FFP
            if transFFP
                if typeof(ffp[c]) <: NLsolve.SolverResults
                    pkt = ffp[c].zero
                else
                    pkt = ffp[:,c]
                end

                coeffsError[j,c] = SphericalHarmonicExpansions.errorTranslation(coeffsError[j,c],pkt)
            end

            # Set coefficients = 0, if they are smaller than the maximum error
            for l = 0:coeffs[1,1].L
                for m = -l:l
                    coeffs[j,c][l,m] = (abs(coeffs[j,c][l,m]) <= coeffsError[j,c][l,m]) ? 0.0 : coeffs[j,c][l,m]
                end
            end

            # New expansion/function
            expansion[j,c] = sphericalHarmonicsExpansion(coeffs[j,c],x,y,z)+0*x;
            func[j,c] = @fastfunc expansion[j,c]
        end
    end

    return coeffs, coeffsError, expansion, func
end
