module matrix
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
   REAL FUNCTION FindDet(matrix, n)
       IMPLICIT NONE
       REAL(kind=r8), DIMENSION(n,n) :: matrix
       INTEGER, INTENT(IN) :: n
       REAL :: m, temp
       INTEGER :: i, j, k, l
       LOGICAL :: DetExists = .TRUE.
       l = 1
       !Convert to upper triangular form
       DO k = 1, n-1
           IF (matrix(k,k) == 0) THEN
               DetExists = .FALSE.
               DO i = k+1, n
                   IF (matrix(i,k) /= 0) THEN
                       DO j = 1, n
                           temp = matrix(i,j)
                           matrix(i,j)= matrix(k,j)
                           matrix(k,j) = temp
                       END DO
                       DetExists = .TRUE.
                       l=-l
                       EXIT
                   ENDIF
               END DO
               IF (DetExists .EQV. .FALSE.) THEN
                   FindDet = 0
                   return
               END IF
           ENDIF
           DO j = k+1, n
               m = matrix(j,k)/matrix(k,k)
               DO i = k+1, n
                   matrix(j,i) = matrix(j,i) - m*matrix(k,i)
               END DO
           END DO
       END DO
       
       !Calculate determinant by finding product of diagonal elements
       FindDet = l
       DO i = 1, n
           FindDet = FindDet * matrix(i,i)
       END DO
       
   END FUNCTION FindDet

   REAL FUNCTION logdet(matrix, n)
       IMPLICIT NONE
       REAL(kind=r8), DIMENSION(n,n) :: matrix
       INTEGER, INTENT(IN) :: n
       REAL :: m, temp
       INTEGER :: i, j, k, l
       LOGICAL :: DetExists = .TRUE.
       l = 1
       !Convert to upper triangular form
       DO k = 1, n-1
           IF (matrix(k,k) == 0) THEN
               DetExists = .FALSE.
               DO i = k+1, n
                   IF (matrix(i,k) /= 0) THEN
                       DO j = 1, n
                           temp = matrix(i,j)
                           matrix(i,j)= matrix(k,j)
                           matrix(k,j) = temp
                       END DO
                       DetExists = .TRUE.
                       l=-l
                       EXIT
                   ENDIF
               END DO
               IF (DetExists .EQV. .FALSE.) THEN
                   logdet = 0
                   return
               END IF
           ENDIF
           DO j = k+1, n
               m = matrix(j,k)/matrix(k,k)
               DO i = k+1, n
                   matrix(j,i) = matrix(j,i) - m*matrix(k,i)
               END DO
           END DO
       END DO
       
       !Calculate determinant by finding product of diagonal elements
       logdet = l
       DO i = 1, n
           logdet = logdet + log(matrix(i,i))
       END DO
       
   END FUNCTION logdet
end module matrix
