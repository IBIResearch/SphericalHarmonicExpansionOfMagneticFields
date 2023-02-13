"""
    magneticField(coords::Array{T,2}, field::Union{Array{T,2},Array{T,3}}, R::T, center::Array{T,1}, L::Int64,
                  x::Variable, y::Variable, z::Variable;
                  calcSolid::Bool=true) where T <: Real
*Description:*  Calculation of the spherical harmonic coefficients and expansion based on the measured t-design\\
 \\
*Input:*
- `coords`      - Coordinates of the measured t-design (size = (3,N))
- `field`       - Measured field (size = (J, N [, C])) with J <= 3
- `R`           - radius of the measured sphere
- `center`      - center of the measured sphere
- `L`           - Order up to which the coeffs be calculated
- `x, y, z`     - Cartesian coordinates
**kwargs:**
- `calcSolid`   - Boolean (default: true)\\
    false -> spherical coefficients\\
    true -> solid coefficients
*Output:*
- `coeffs`    - spherical/solid coefficients, type: Array{SphericalHarmonicCoefficients}(3,C)
- `expansion` - related expansion (Cartesian polynomial), type: Array{Polynomial}(3,C)
- `func`      - expansion converted to a function, type: Array{Function}(3,C)
"""
function magneticField(coords::Array{T,2}, field::Union{Array{T,2},Array{T,3}}, R::T, center::Vector{T}, L::Int64,
		       x::Variable, y::Variable, z::Variable;
		       calcSolid::Bool=true) where T <: Real

  # transpose coords if its dimensions do not fit
  if size(coords,1) != 3
    coords = coords'
  end

  # test dimensions of field array
  if size(field,1) > 3
    throw(DimensionMismatch("The measured field has more than 3 entries in the first dimension: $(size(field,1))"))
  elseif size(field,2) != size(coords,2)
    throw(DimensionMismatch("The field vector does not match the size of the tdesign: $(size(field,2)) != $(size(coords,2))"))
  end

  func= Array{Function}(undef,size(field,1),size(field,3))
  expansion = Array{Polynomial}(undef,size(field,1),size(field,3))
  coeffs = Array{SphericalHarmonicCoefficients}(undef,size(field,1),size(field,3))

  # rescale coordinates to t-design on unit sphere
  coords = coords .- center
  coords *= 1/R
  for c = 1:size(field,3)
    # Calculation of the coefficients
    for j = 1:size(field,1)

        coeffs[j,c] = SphericalHarmonicExpansions.sphericalQuadrature(field[j,:,c],coords',L);
        coeffs[j,c].R = R
	
        normalize!(coeffs[j,c],R)
	
	# convert spherical into solid coefficients
        if calcSolid
            solid!(coeffs[j,c])
        end

        # calculation of the expansion
        expansion[j,c] = sphericalHarmonicsExpansion(coeffs[j,c],x,y,z) + 0*x;
        func[j,c] = @fastfunc expansion[j,c]+0*x+0*y+0*z
    end
  end

  return coeffs, expansion, func
end
