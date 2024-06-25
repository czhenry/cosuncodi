export permutation, B, B_hat, G, G_hat, J, J_hat, M, M_hat

function permutation(N, M)
    factorial(N)//factorial(N-M)
end

function B(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in 0:N; for j in 0:N
      if (i<j)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))

        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function B_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in 0:N; for j in 0:N
      if (i<j)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^j*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function G(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function G_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-N-1)*binomial(big(i),big(i- N-1))*permutation(big(j), big(i-N-1))*factorial(big(N+1))
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-N-1)*binomial(i,i- N-1)*permutation(j, i-N-1)*factorial(N+1)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(i-N-1)*binomial(i, i-N-1)*permutation(j, i-N-1)*factorial(N+1)
      end
    end
    end
    result
end

function J(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(N+1)*binomial(big(i), big(N+1))*permutation(big(j), big(i-N-1))*factorial(big(N+1))
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(N+1)*binomial(i, N+1)*permutation(j, i-N-1)*factorial(N+1)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(N+1)*binomial(i, N+1)*permutation(j, i-N-1)*factorial(N+1)
      end
    end
    end
    result
end

function J_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^j*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function M(numType, N::Int64, a)
    result=Matrix{numType}(undef, 2*N+2, N+1)
    for i in 0:(2*N+1); for j in 0:N
      if (i>(j+N+1))||(j>i)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(N+1), big(i-j))
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(N+1, i-j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-j)*binomial(N+1, i-j)
      end
    end
    end
    result
end

function M_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, 2*N+2, N+1)
    for i in 0:(2*N+1); for j in 0:N
      if (i>(j+N+1))||(i<(N+1))
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-N-1)*binomial(big(j), big(i-N-1))
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-N-1)*binomial(j, i-N-1)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-N-1)*binomial(j, i-N-1)
      end
    end
    end
    result
end
