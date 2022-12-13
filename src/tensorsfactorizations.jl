#= # This script has been copy-pasted from https://github.com/mhauru/TensorFactorizations.jl

"""
    tensorsplit(A, a, b; kwargs...)

Calls tensorsvd with the arguments given to it to decompose the given tensor
A with indices a on one side and indices b on the other.  It then splits
the diagonal matrix of singular values into two with a square root and
multiplies these weights into the isometric tensors.  Thus tensorsplit ends
up splitting A into two parts, which are then returned, possibly together
with auxiliary data such as a truncation error. If the keyword argument
hermitian=true, an eigenvalue decomposition is used in stead of an SVD. All
the keyword arguments are passed to either tensorsvd or tensoreig.

See tensorsvd and tensoreig for further documentation.
"""
function tensorsplit(args...; kwargs...)
    # Find the keyword argument hermitian.
    # TODO This is awful, why do I have to do this?
    hermitian = false
    for (key, value) in kwargs
        key == :hermitian && (hermitian = value)
    end

    if hermitian
        res = tensoreig(args...; kwargs...)
        S, U = res[1:2]
        Vt_perm = [ndims(U), (1:(ndims(U) - 1))...]
        Vt = conj!(tensorcopy(U, collect(1:ndims(U)), Vt_perm))
        S = Diagonal(S)
        if !isposdef(S)
            S = complex.(S)
        end
        auxdata = res[3:end]
    else
        res = tensorsvd(args...; kwargs...)
        U, S, Vt = res[1:3]
        S = Diagonal(S)
        auxdata = res[4:end]
    end
    S_sqrt = sqrt.(S)
    A1 = tensorcontract(U, (1:(ndims(U) - 1)..., :a), S_sqrt, (:a, :b))
    A2 = tensorcontract(S_sqrt, (:b, :a), Vt, (:a, 1:(ndims(Vt) - 1)...))
    return A1, A2, auxdata...
end =#

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

#= """
    tensorsvd(A, a, b;
              chis=nothing, eps=0,
              return_error=false, print_error=false,
              break_degenerate=false, degeneracy_eps=1e-6,
              norm_type=:frobenius)

Singular valued decomposes a tensor A. The indices of A are
permuted so that the indices listed in the Array/Tuple a are on the "left"
side and indices listed in b are on the "right".  The resulting tensor is
then reshaped to a matrix, and this matrix is SVDed into U*diagm(S)*Vt.
Finally, the unitary matrices U and Vt are reshaped to tensors so that
they have a new index coming from the SVD, for U as the last index and for
Vt as the first, and U has indices a as its first indices and V has
indices b as its last indices.

If eps>0 then the SVD may be truncated if the relative error can be kept
below eps. For this purpose different dimensions to truncate to can be tried,
and these dimensions should be listed in chis. If chis is nothing (the
default) then the full range of possible dimensions is tried. If
break_degenerate=false (the default) then the truncation never cuts between
degenerate singular values. degeneracy_eps controls how close the values need
to be to be considered degenerate.

norm_type specifies the norm used to measure the error. This defaults to
:frobenius, which means that the error measured is the Frobenius norm of the
difference between A and the decomposition, divided by the Frobenius norm of
A.  This is the same thing as the 2-norm of the singular values that are
truncated out, divided by the 2-norm of all the singular values. The other
option is :trace, in which case a 1-norm is used instead.

If print_error=true the truncation error is printed. The default is false.

If return_error=true then the truncation error is also returned.
The default is false.

Note that no iterative techniques are used, which means choosing to truncate
provides no performance benefits: The full SVD is computed in any case.

Output is U, S, Vt, and possibly error. Here S is a vector of
singular values and U and Vt are isometric tensors (unitary if the matrix
that is SVDed is square and there is no truncation) such that  U*diag(S)*Vt =
A, up to truncation errors.
"""
function tensorsvd(
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
)
    # Create the matrix and SVD it.
    A, shp_a, shp_b = to_matrix(A, a, b; return_tensor_shape=true)
    fact = svd(A)
    U, S, Vt = fact.U, fact.S, fact.Vt

    # Find the dimensions to truncate to and the error caused in doing so.
    chi, error = find_trunc_dim(S, chis, eps, break_degenerate, degeneracy_eps, norm_type)
    # Truncate
    S = S[1:chi]
    U = U[:, 1:chi]
    Vt = Vt[1:chi, :]

    if print_error
        println("Relative truncation error ($norm_type norm) in SVD: $error")
    end

    # Reshape U and V to tensors with shapes matching the shape of A and
    # return.
    dim = size(S)[1]
    U_tens = reshape(U, shp_a..., dim)
    Vt_tens = reshape(Vt, dim, shp_b...)
    retval = (U_tens, S, Vt_tens)
    if return_error
        retval = (retval..., error)
    end
    return retval
end
 =#
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
    A = tensorcopy(A, collect(1:ndims(A)), perm)
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
