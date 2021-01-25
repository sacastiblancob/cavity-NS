!               ***********************
                MODULE DECLARATIONS_CSC
!               ***********************
!
!
!***********************************************************************
! STRUCTURE DECLARATION FOR CSC STRUCTURES
!***********************************************************************
!
!brief    STRUCTURE CSC FOR STORAGE VALUES, ROWS INDICES AND COLUMN
!         STARTS
!
!history  Sergio Castiblanco
!+        16/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE CSC_OBJ
!
!  NUMBER OF NONZEROS
       INTEGER :: NZ
!
!  NUMBER OF COLUMNS
       INTEGER :: NC
!
!  NUMBER OF ROWS
       INTEGER :: NR
!
!  NAME
       CHARACTER(LEN=6) :: NAM
!
!  VALUES
       DOUBLE PRECISION, POINTER, DIMENSION(:) :: V => NULL()
!
!  ROW INDICES
       INTEGER, POINTER, DIMENSION(:) :: R => NULL()
!
!  COLUMN RELATED INDICES
!  (WHAT IS THE POSITION IN "V" WHERE THE COLUMN START)
       INTEGER, POINTER, DIMENSION(:) :: C => NULL()
!
      END TYPE CSC_OBJ
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END MODULE DECLARATIONS_CSC



