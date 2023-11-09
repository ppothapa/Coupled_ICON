MODULE mo_ocean_math_utils

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: solve_tridiag_block

  INTERFACE solve_tridiag_block
    MODULE PROCEDURE solve_tridiag_block_wp
  END INTERFACE

CONTAINS

!>
!! Solve a set of tridiagonal systems.
!! Each system is formed by `a(i,:)`, `b(i,:)`, `c(i,:)`, `d(i,:)`. The number of rows in each
!! system is given by `n(i)`:
!! \f[ A_{i,k}x_{i,k-1} + B_{i,k}x_{i,k} + C_{i,k}x{i,k+1} = D_{i,k}, \quad 1 \le k \le n_i \f]
!!
!! The routine performs Gaussian elimination on the systems, followed by backsubstitution to
!! compute the solution. The direction of elimination can be chosen via `eliminate_upper` as
!! first to last row (`eliminate_upper=.false.`) or last to first row.
!!
!! \notice Although they are irrelevant to the solution, `a(i,1)` and `c(i,n(i))` have to be finite
!! floating-point values.
SUBROUTINE solve_tridiag_block_wp (a, b, c, d, x, n, eliminate_upper, use_acc)
    REAL(wp), INTENT(IN) :: a(:,:) !< Subdiagonal entries (nsys,nrow).
    REAL(wp), INTENT(IN) :: b(:,:) !< Diagonal entries (nsys,nrow).
    REAL(wp), INTENT(IN) :: c(:,:) !< Superdiagonal entries (nsys,nrow).
    REAL(wp), INTENT(IN) :: d(:,:) !< Right-hand side (nsys,nrow).
    REAL(wp), INTENT(INOUT) :: x(:,:) !< Solutions (nsys,nrow).
    INTEGER, INTENT(IN) :: n(:) !< Number of rows (nsys).
    LOGICAL, INTENT(IN) :: eliminate_upper !< Solve by eliminating the upper diagonal, not the lower.
    LOGICAL, OPTIONAL, INTENT(IN) :: use_acc

    REAL(wp) :: cp(SIZE(a,1), SIZE(a,2))
    REAL(wp) :: dp(SIZE(a,1), SIZE(a,2))
    REAL(wp) :: fxa
    integer :: i, j, maxn
    LOGICAL :: lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !$ACC DATA CREATE(cp, dp) IF(lacc)

    maxn = MAXVAL(n)

    IF (eliminate_upper) THEN

      ! initialize a-prime (cp) and d-prime
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO j = 1, SIZE(a,1)
        IF (n(j) > 0) THEN
          cp(j,n(j)) = a(j,n(j)) / b(j,n(j))
          dp(j,n(j)) = d(j,n(j)) / b(j,n(j))
        END IF
      END DO

      ! solve for vectors a-prime and d-prime
      !$ACC LOOP SEQ
      DO i = maxn-1,1,-1
        !$ACC LOOP GANG VECTOR
        DO j = 1, SIZE(a,1)
          IF (i <= n(j)-1) THEN
            fxa = 1.0_wp / (b(j,i) - cp(j,i+1) * c(j,i))
            cp(j,i) = a(j,i) * fxa
            dp(j,i) = (d(j,i) - dp(j,i+1) * c(j,i)) * fxa
          END IF
        END DO
      END DO

      ! initialize x
      !$ACC LOOP GANG VECTOR
      DO j = 1, SIZE(a,1)
        IF (n(j) > 0) THEN
          x(j,1) = dp(j,1)
        END IF
      END DO

      ! solve for x from the vectors a-prime and d-prime
      !$ACC LOOP SEQ
      DO i = 2, maxn
        !$ACC LOOP GANG VECTOR
        DO j = 1, SIZE(a,1)
          IF (i <= n(j)) THEN
            x(j,i) = dp(j,i) - cp(j,i)*x(j,i-1)
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)

    ELSE ! eliminate_upper

      ! initialize c-prime and d-prime
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO j = 1, SIZE(a,1)
        IF (n(j) > 0) THEN
          cp(j,1) = c(j,1) / b(j,1)
          dp(j,1) = d(j,1) / b(j,1)
        END IF
      END DO

      ! solve for vectors c-prime and d-prime
      !$ACC LOOP SEQ
      DO i = 2, maxn
        !$ACC LOOP GANG VECTOR
        DO j = 1, SIZE(a,1)
          IF (i <= n(j)) THEN
            fxa = 1.0_wp / (b(j,i) - cp(j,i-1) * a(j,i))
            cp(j,i) = c(j,i) * fxa
            dp(j,i) = (d(j,i) - dp(j,i-1) * a(j,i)) * fxa
          END IF
        END DO
      END DO

      ! initialize x
      !$ACC LOOP GANG VECTOR
      DO j = 1, SIZE(a,1)
        IF (n(j) > 0) THEN
          x(j,n(j)) = dp(j,n(j))
        END IF
      END DO

      ! solve for x from the vectors c-prime and d-prime
      !$ACC LOOP SEQ
      DO i = maxn-1, 1, -1
        !$ACC LOOP GANG VECTOR
        DO j = 1, SIZE(a,1)
          IF (i <= n(j) - 1) THEN
            x(j,i) = dp(j,i) - cp(j,i)*x(j,i+1)
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA
  END SUBROUTINE solve_tridiag_block_wp

END MODULE mo_ocean_math_utils 
