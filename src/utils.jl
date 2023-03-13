function correlation(a::Array, b=nothing::Union{Nothing,Array}; trunc=nothing::Union{Nothing,Int})
    fra = fft(a)
    if isnothing(b)
        sf = conj.(fra) .* fra
    else
        if length(b)!= length(a)
            throw(ArgumentError("Argument should be of the same size"))
        end
        # frb = fft(vcat(b,zeros(length(a)-length(b))))
        frb = fft(b)
        sf = conj.(frb) .* fra
    end

    if isnothing(trunc)
        len_trunc = length(a)
    else
        len_trunc = min(length(a), trunc+1)
    end
    res = ifft(sf)
    cor = real.(res[1:len_trunc]) ./ range(start=length(a),stop=length(a) - len_trunc+1,step=-1 )
    return cor
end
