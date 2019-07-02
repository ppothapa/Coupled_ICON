INTERFACE
SUBROUTINE RRTM_TAUMOL12 (KIDIA,KFDIA,KLEV,P_TAU,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,P_ONEMINUS,&
 & P_COLH2O,P_COLCO2,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC,&
 & PRAT_H2OCO2, PRAT_H2OCO2_1) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE PARRRTM , ONLY : JPBAND
USE YOERRTM , ONLY : JPGPT ,NG12 ,NGS11
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB) ,INTENT(OUT) :: P_TAU(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC11(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JP(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_ONEMINUS
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLH2O(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLCO2(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_LAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_INDSELF(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PFRAC(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PRAT_H2OCO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PRAT_H2OCO2_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN) :: K_INDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FORFRAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FORFAC(KIDIA:KFDIA,KLEV)
END SUBROUTINE RRTM_TAUMOL12
END INTERFACE
