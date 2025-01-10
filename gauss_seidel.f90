program gauss_seidel
    implicit none   ! Sem declarações implícitas

    ! Variáveis
    ! Aqui, partimos da definição Ax = b para sistemas lineares, em que a é uma matriz n x n, x é o chute inicial n x 1 e b é uma matriz de constantes n x 1.
    integer :: n, i, j, k, max_iterations
    real :: tolerance, error
    real, allocatable :: A(:,:), b(:), x(:), x_old(:)   ! Vetores alocados na heap, tamanho determinado no runtime (necessário defini-los assim em Fortran)

    ! Tamanho do sistema linear (obrigatoriamente se torna uma matriz quadrada, deve seguir a definição de nossa representação)
    print *, "Enter the column size of the coefficient matrix (n):"
    read *, n

    ! Valida se o número inserido é maior ou igual a 1
    if (n < 1) then
        print *, "Matrix size must be at least 1."
        stop
    end if

    ! Aloca a memória correspondente ao tamanho inserido para as matrizes
    allocate(A(n, n))
    allocate(b(n))
    allocate(x(n))
    allocate(x_old(n))

    ! Inputs da matriz de coeficientes A
    print *, "Enter the coefficients of matrix A row by row:"
    do i = 1, n
        do j = 1, n
            write(*, '(A, I2, A, I2, A)') "A(", i, ",", j, "):"
            read *, A(i, j)
        end do
    end do

    ! Inputs do vetor b de constantes
    print *, "Enter the constants vector b:"
    do i = 1, n
        write(*, '(A, I2, A)') "b(", i, "):"
        read *, b(i)
    end do

    ! Input do vetor x(0) (chute inicial)
    print *, "Enter the initial guess vector x:"
    do i = 1, n
        write(*, '(A, I2, A)') "x(", i, "):"
        read *, x(i)
    end do

    ! Tolerância e máximo de iterações (para não estourar a call-stack)
    print *, "Enter the tolerance for convergence:"
    read *, tolerance

    print *, "Enter the maximum number of iterations:"
    read *, max_iterations

    error = 1.0  ! Inicializa a margem de erro para um valor grande arbitrário (1.0)
    k = 0        ! Conta as iterações (k)

    ! Algoritmo de Gauss-Seidel, executado enquanto o erro é maior do que a tolerância e k é menor do que o máximo de iterações
    do while (error > tolerance .and. k < max_iterations)
        k = k + 1          ! Incrementa contador de iterações
        x_old = x          ! Armazena os valores de x em x_old (valores de x(k + 1) em x(k) para a próxima iteração (aqui entra o aspecto recursivo)

        ! Calcula cada termo do novo vetor solução (x) seguindo a separação pela diagonal de A
        do i = 1, n
            x(i) = b(i)     ! Atribuímos b(i) a x(i) para inicializar a variável
            do j = 1, n
                if (j /= i) then    ! Se não for um elemento da diagonal
                    x(i) = x(i) - A(i, j) * x(j)  ! Aplica x(k + 1) = b(i) - A(i, j) * x(k + 1)(j)
                end if
            end do
            x(i) = x(i) / A(i, i)  ! Divide o termo calculado pela diagonal
        end do

        ! Calcula o erro absoluto (maior valor absoluto das diferenças entre os termos correspondentes de cada vetor)
        error = maxval(abs(x - x_old))

        ! Debug de cada iteração
        print *, "Iteration", k, ": x =", x, ", Error =", error
    end do

    ! Se  o erro absoluto for menor ou igual do que a tolerância, converge
    if (error <= tolerance) then
        print *, "Gauss-Seidel method converged in", k, "iterations."
        print *, "Solution vector x:"
        write(*, '(100(F10.5))') x
    else
        print *, "Gauss-Seidel method did not converge within the maximum iterations."
    end if

    ! Libera espaço alocado na heap
    deallocate(A)
    deallocate(b)
    deallocate(x)
    deallocate(x_old)

end program gauss_seidel