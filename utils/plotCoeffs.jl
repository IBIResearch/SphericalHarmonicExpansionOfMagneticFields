"""
    fignumber = plotCoeffs(coeffs::Array{SphericalHarmonicCoefficients}, L::Int=2, c::Int=1;
                    R::T=1.0, fignumber::Int=1,
                    separate::Bool=false, sepnorm::Bool=true,
                    ymin::T=0.0, ymax::T=0.0, fontsize::Real=12,
                    coeffsError=0, plotError::Bool=false,
                    useMilli::Bool=false, normalizemm::Bool=false) where T <: Real
*Description:*  Plot of the coefficients\\
 \\
*Input:*
- `coeffs`      - Coefficients
- `c`           - Index (of the current) of the coeffs, which should be plotted (default: 1)
- `L`           - Order of the largest plotted coeffs (default: 2)
**kwargs:**
- `R`           - Normalization radius (default: 1.0)\\
                  -> if != 1.0: normalize coeffs with 1/R
- `fignumber`   - Number of the figure   (default: 1)
- `separate`    - Boolean -> only for length(c)=1 (default: false)\\
                  true -> 3 subplots for x,y,z will be generated\\
                  false -> 1 plot with all coeffs
- `sepnorm`	- Boolean -> only for seperate=true and ymin=ymax (default: true)\\
		  true -> normalized y-axis to same min/max for all 3 subplots 
- `ymin, ymax`  - Minimum and maximum for yticks (default: 0)\\
                    (if ymin == ymax yticks are not set manually)
- `fontsize`    - Fontsize (default: 12)
- `coeffsError` - Error of the coeffs (default: 0)
- `plotError`   - Boolean (default: false)\\
                  true -> error of the coefficients will be plotted (-> coeffsError)
- `useMilli`    - Boolean (default: false)\\
                  true -> convert tesla to millitesla
- `normalizemm` - Boolean (default: false)\\
                  true -> normalize to R = 1 mm
*Output:*
- `fignumber` = fignumber + 1
"""
function plotCoeffs(coeffs::Array{SphericalHarmonicCoefficients}, L::Int=2, c::Int=1;
                    R::T=1.0, fignumber::Int=1,
                    separate::Bool=false, sepnorm::Bool=true,
		    ymin::T=0.0, ymax::T=0.0, fontsize::Real=12,
                    coeffsError=0, plotError::Bool=false,
                    useMilli::Bool=false, normalizemm::Bool=false) where T <: Real

  # Labels of x-axis
  xx = Array{String}(undef,(coeffs[1,c[1]].L+1)^2)
  i = 1
  for l=0:coeffs[1,c[1]].L
    for m=-l:l
      xx[i] = "[$l,$m]"
      i += 1
    end
  end
  ii = 1:1:(coeffs[1,c[1]].L+1)^2

  # Field dimension
  J = size(coeffs,1)

  # convert field values of the plot to mT if useMilli = true
  convMilli = useMilli ? 1000 : 1;

  # normalize coefficients to R = 0.001 (= 1 mm) if normalizemm = true
  coeffsR = deepcopy(coeffs)
  coeffsErrorR = deepcopy(coeffsError)

  coeffsR = normalizemm ? normalize.(coeffsR,coeffsR[1,c[1]].R*1000) : coeffsR
  coeffsErrorR = (normalizemm && plotError) ? normalize.(coeffsErrorR,coeffsErrorR[1,c[1]].R*1000) : coeffsErrorR

  # if R != 1.0 normalize coefficients
  if R != 1.0
    normalize!.(coeffsR,1/R)
    if plotError
      normalize!.(coeffsErrorR,1/R)
    end
  end

  # set text for ylabel
  if coeffsR[1,c[1]].R == 1 && coeffsR[1,c[1]].solid && useMilli
    ytext = L"$\gamma_{l,m}\;/\;mT/m^l$"
  elseif coeffsR[1,c[1]].R == 1 && coeffsR[1,c[1]].solid && !useMilli
    ytext = L"$\gamma_{l,m}\;/\;T/m^l$"
  elseif coeffsR[1,c[1]].R == 1 && !(coeffsR[1,c[1]].solid) && useMilli
    ytext = L"$c_{l,m}\;/\;mT/m^l$"
  elseif coeffsR[1,c[1]].R == 1 && !(coeffsR[1,c[1]].solid) && !useMilli
    ytext = L"$c_{l,m}\;/\;T/m^l$"
  elseif coeffsR[1,c[1]].R == 0.001 && coeffsR[1,c[1]].solid && useMilli
    ytext = L"$\gamma_{l,m}\;/\;mT/mm^l$"
  elseif coeffsR[1,c[1]].R == 0.001 && coeffsR[1,c[1]].solid && !useMilli
    ytext = L"$\gamma_{l,m}\;/\;T/mm^l$"
  elseif coeffsR[1,c[1]].R == 0.001 && !(coeffsR[1,c[1]].solid) && useMilli
    ytext = L"$c_{l,m}\;/\;mT/mm^l$"
  elseif coeffsR[1,c[1]].R == 0.001 && !(coeffsR[1,c].solid) && !useMilli
    ytext = L"$c_{l,m}\;/\;T/mm^l$"
  elseif coeffsR[1,c[1]].R != 1 && coeffsR[1,c[1]].solid && useMilli
    ytext = L"$\gamma_{l,m}^R\;/\;mT$"
  elseif coeffsR[1,c[1]].R != 1 && coeffsR[1,c[1]].solid && !useMilli
    ytext = L"$\gamma_{l,m}^R\;/\;T$"
  elseif coeffsR[1,c[1]].R != 1 && !(coeffsR[1,c[1]].solid) && useMilli
    ytext = L"$c_{l,m}^R\;/\;mT$"
  else # coeffsR[1,c].R != 1 && !(coeffsR[1,c].solid) && !useMilli
    ytext = L"$c_{l,m}^R\;/\;T$"
  end
	
  # set color vector (use same colors for separate and combined plots)
  # blue, green, yellow
  colorvec = [(0.0,0.2862,0.5725,1.0), (0.5412,0.7412,0.1412,1.0), (1.0,0.8745,0.0,1.0)] 
	

  LL = (L+1)^2
  if separate
        ## Plot in 3 subplots ##
        figure(fignumber,figsize=(14,15))
        ti = ["x","y","z"]
	
	# set consistent min and max for y-axis
	if sepnorm && ymin == ymax
	  ymin = minimum([minimum(coeffsR[j,c].c[1:LL].*convMilli) for j=1:3])
	  ymax = maximum([maximum(coeffsR[j,c].c[1:LL].*convMilli) for j=1:3])
	  ydiff = maximum([abs(ymin),abs(ymax)])
	  ymin -= 0.1*ydiff
	  ymax += 0.1*ydiff
	end

        for j = 1:J
            subplot(J,1,j)
            bar(ii[1:LL],(coeffsR[j,c].c[1:LL]).*convMilli,color=colorvec[j],label=ti[j])
            if plotError
		# add error lines
                plot(ii[1:LL],abs.(coeffsErrorR[j,c].c[1:LL]).*convMilli,color=colorvec[j],alpha=0.5,linewidth="2",label="Error")
                plot(ii[1:LL],-abs.(coeffsErrorR[j,c].c[1:LL]).*convMilli,color=colorvec[j],alpha=0.5,linewidth="2")
            end

	    # labels, ticks, ...
            xlabel("[l,m]",fontsize=fontsize+2)
            ylabel(ytext,fontsize=fontsize+2)

            grid(true)

            xticks(ii[1:LL],xx[1:LL],fontsize=fontsize)
            yticks(fontsize=fontsize)
            if ymin != ymax
                ylim(ymin,ymax)
            end
            legend(fontsize=fontsize)

            title(ti[j],fontsize=fontsize)
        end
        tight_layout()
  else
        ## Plot of the coefficients in a single plot ##
        figure(fignumber,figsize=(5*L,5))

        bar(ii[1:LL].-0.3,(coeffsR[1,c].c[1:LL]).*convMilli,color=colorvec[1],label="x",width=0.3)
        J >= 2 && bar(ii[1:LL],(coeffsR[2,c].c[1:LL]).*convMilli,color=colorvec[2],label="y",width=0.3)
        J == 3 && bar(ii[1:LL].+0.3,(coeffsR[3,c].c[1:LL]).*convMilli,color=colorvec[3],label="z",width=0.3)
        if plotError
	    # add error lines
            plot(ii[1:LL].-0.3,abs.(coeffsErrorR[1,c].c[1:LL]).*convMilli,color=colorvec[1],alpha=0.5,linewidth="2",label="Error in x")
            J >= 2 && plot(ii[1:LL],abs.(coeffsErrorR[2,c].c[1:LL]).*convMilli,color=colorvec[2],alpha=0.5,linewidth="2",label="Error in y")
            J == 3 && plot(ii[1:LL].+0.3,abs.(coeffsErrorR[3,c].c[1:LL]).*convMilli,color=colorvec[3],alpha=0.5,linewidth="2",label="Error in z")
            plot(ii[1:LL].-0.3,-abs.(coeffsErrorR[1,c].c[1:LL]).*convMilli,color=colorvec[1],alpha=0.5,linewidth="2")
            J >= 2 && plot(ii[1:LL],-abs.(coeffsErrorR[2,c].c[1:LL]).*convMilli,color=colorvec[2],alpha=0.5,linewidth="2")
            J == 3 && plot(ii[1:LL].+0.3,-abs.(coeffsErrorR[3,c].c[1:LL]).*convMilli,color=colorvec[3],alpha=0.5,linewidth="2")
        end

	# labels, ticks, ...
        xlabel("[l,m]",fontsize=fontsize+2)
        ylabel(ytext,fontsize=fontsize+2)

        grid(true)

        xticks(ii[1:LL],xx[1:LL],fontsize=fontsize)
        yticks(fontsize=fontsize)
        if ymin != ymax
            ylim(ymin,ymax)
            # yticks(linspace(ymin,ymax,5))
        end
        legend(fontsize=fontsize)

        tight_layout()
  end

  return fignumber+1
end
