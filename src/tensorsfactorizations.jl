# This script has been copy-pasted from https://github.com/mhauru/TensorFactorizations.jl

"""
    tensoreig(A, a, b; chis=nothing, eps=0,
              return_error=false, print_error=false,
              break_degenerate=false, degeneracy_eps=1e-6,
              norm_type=:frobenius, hermitian=false)

Finds the "right" eigenvectors and eigenvalues of A. The indices of A are
permuted so that the indices listed in the Array/Tuple a are on the "left"
side and indices listed in b are on the "right".  The resulting tensor is
then reshaped to a matrix, and eig is called on this matrix to get the vector
of eigenvalues E and matrix of eigenvectors U. Finally, U is reshaped to
a tensor that has as its last index the one that enumerates the eigenvectors
and the indices in a as its first indices.

Truncation and error printing work as with tensorsvd.

Note that no iterative techniques are used, which means that choosing to
truncate provides no performance benefits: All the eigenvalues are computed
in any case.

The keyword argument hermitian (false by default) tells the algorithm
whether the reshaped matrix is Hermitian or not. If hermitian=true, then A =
U*diagm(E)*U' up to the truncation error.

Output is E, U, and possibly error, if return_error=true. Here E is a
vector of eigenvalues values and U[:,...,:,k] is the kth eigenvector.
"""
function tensoreig(
    A,
    a,
    b;
    chis=nothing,
    eps=0,
    return_error=false,
    print_error=false,
    break_degenerate=false,
    degeneracy_eps=1e-6,
    norm_type=:frobenius,
    hermitian=false,
)
    # Create the matrix and decompose it.
    A, shp_a, shp_b = to_matrix(A, a, b; return_tensor_shape=true)
    if hermitian
        A = Hermitian((A + A') / 2)
    end
    fact = eigen(A)
    E, U = fact.values, fact.vectors
    # Sort the by largest magnitude eigenvalue first.
    perm = sortperm(abs.(E); rev=true)
    if perm != collect(1:length(E))
        E = E[perm]
        U = U[:, perm]
    end

    # Find the dimensions to truncate to and the error caused in doing so.
    chi, error = find_trunc_dim(E, chis, eps, break_degenerate, degeneracy_eps, norm_type)

    # Truncate
    E = E[1:chi]
    U = U[:, 1:chi]

    if print_error
        println("Relative truncation error ($norm_type norm) in eig: $error")
    end

    # Reshape U to a tensor with a shape matching the shape of A and
    # return.
    dim = size(E)[1]
    U_tens = reshape(U, shp_a..., dim)
    retval = (E, U_tens)
    if return_error
        retval = (retval..., error)
    end
    return retval
end

"""
Format the bond dimensions listed in chis to a standard format.
"""
function format_trunc_chis(v, chis, eps)
    max_dim = length(v)
    if chis == nothing
        if eps > 0
            # Try all possible chis.
            chis = collect(1:max_dim)
        else
            # No truncation.
            chis = [max_dim]
        end
    else
        if isa(chis, Number)
            # Wrap an individual number in an Array.
            chis = [chis]
        else
            # Make sure chis is an Array, and not, say, a tuple.
            chis = collect(chis)
        end
        if eps == 0
            chis = [maximum(chis)]
        else
            sort!(chis)
        end
    end
    # If some of the chis are larger than max_dim, get rid of them.
    for (i, chi) in enumerate(chis)
        if chi >= max_dim
            chis[i] = max_dim
            chis = chis[1:i]
            break
        end
    end
    return chis
end

"""
Transpose A so that the indices listed in a are on the left and the indices
listed in b on the right, and reshape A into a matrix.

a and b should be arrays of Integers that together include the numbers from 1
to ndims(A). Alternatively they can be just individual Integers, if only one
index is left on one side.

If return_tensor_shape is true (by default it's not) return, in addition to the
matrix, the shape of the tensor after the transpose but before the reshape.
"""
function to_matrix(A, a, b; return_tensor_shape=false)
    # Make sure a and b are Arrays.
    if isa(a, Number)
        a = [a]
    elseif !isa(a, Array)
        a = collect(a)
    end
    if isa(b, Number)
        b = [b]
    elseif !isa(b, Array)
        b = collect(b)
    end

    # Permute the indices of A to the right order
    perm = vcat(a, b)
    A = tensorcopy(collect(1:ndims(A)), A, perm)
    # The lists shp_a and shp_b list the dimensions of the bonds in a and b
    shp = size(A)
    shp_a = shp[1:length(a)]
    shp_b = shp[(end + 1 - length(b)):end]

    # Compute the dimensions of the the matrix that will be formed when
    # indices of a and b are joined together.
    dim_a = prod(shp_a)
    dim_b = prod(shp_b)

    A = reshape(A, dim_a, dim_b)

    if return_tensor_shape
        return A, shp_a, shp_b
    else
        return A
    end
end

"""
Finds the dimension to which v should be truncated, and the error caused in
this truncation. See documentation for tensorsvd for the meaning of the
different arguments.
"""
function find_trunc_dim(
    v, chis, eps, break_degenerate=false, degeneracy_eps=1e-6, norm_type=:frobenius
)
    # Put chis in a standard format.
    chis = format_trunc_chis(v, chis, eps)
    v = abs.(v)
    if !issorted(v; rev=true)
        sort!(v; rev=true)
    end

    if norm_type == :frobenius
        v = v .^ 2
        eps = eps^2
    elseif norm_type == :trace
        # Nothing to be done
    else
        throw(ArgumentError("Unknown norm_type $norm_type."))
    end
    sum_all = sum(v)
    # Find the smallest chi for which the error is small enough.
    # If none is found, use the largest chi.
    if sum_all != 0
        for i in 1:length(chis)
            chi = chis[i]
            if !break_degenerate
                # Make sure that we don't break degenerate singular values by
                # including one but not the other by decreasing chi if
                # necessary.
                while 0 < chi < length(v)
                    last_in = v[chi]
                    last_out = v[chi + 1]
                    rel_diff = abs(last_in - last_out) / last_in
                    if rel_diff < degeneracy_eps
                        chi -= 1
                    else
                        break
                    end
                end
            end
            sum_disc = sum(v[(chi + 1):end])
            error = sum_disc / sum_all
            if error <= eps
                break
            end
        end
        if norm_type == :frobenius
            error = sqrt(error)
        end
    else
        error = 0
        chi = minimum(chis)
    end
    return chi, error
end
