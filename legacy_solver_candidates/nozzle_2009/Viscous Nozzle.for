! ___________________________________________________________________
!|                                                                    |
!|  program Nozzle.FOR                                                |
!|                                                                    |
!|         : written for micro-Poissellet flows by DSMC method        |
!|         : 2D version 1.0                                           |
!|         : August, 2008                                             |
!|                                                                    |
!|									                                |
!|____________________________________________________________________|
!======================================================================
      PROGRAM DSMC2
      
	Include 'common.txt'
	include 'property.txt' 
	
	call startcode    
	
	!CALL PLOT_GRID_2_ZONE  

      do while (NPR<NPT)

	NPR=NPR+1

      IF (NPR.LE.NPS) CALL SAMPI2     

      DO J=1,NSP

       DO I=1,NIS
        TIME=TIME+DTM

     
        CALL MOVE2
        CALL INDEXM
        CALL COLLMR
	 end do

        CALL SAMPLE2

	end do
		
	WRITE (*,*) NPR,NM,IDINT(NCOL),TIME

	IPLOT=IPLOT+1
	open(3,FILE='HISTORY.PLT',position="append")
	call cal_var
	
	write(3,'(I8,9F14.8)')IPLOT,Rm0,RmL,u0,ul
	close(3)	

      IF((NPR-(NPR/10)*10)==0) call writerestart
	
	IF((NPR-(NPR/10)*10)==0) CALL OUT2
	
	IF((NPR-(NPR/NSQ)*NSQ)==0) CALL WallHeatFlux

	end do
      
      END program
!========================================================================
	SUBROUTINE startcode

	include 'common.txt'
	
	WRITE (*,*) ' INPUT 0,1 FOR CONTINUING,NEW CALCULATION:- '
      READ (*,*) NQL
!      WRITE (*,*) ' INPUT 0,1 FOR CONTINUING,NEW SAMPLE:- '
!      READ (*,*) NQLS
	
	IF (NQL.EQ.0) THEN
	NQLS=0
	ELSE 
	NQLS=1
	END IF

	IPLOT=0

	CALL SYSTEM ("del *.obj ")
      CALL SYSTEM ("del *.plg ")
	CALL SYSTEM ("del *.opt ")
	CALL SYSTEM ("del *.dsp ")
	CALL SYSTEM ("del *.dsw ")

	IF (NQL.EQ.1) THEN

      CALL INIT2

	open(3,FILE='HISTORY.PLT')
	write(3,*)'Variables= "Print Cycle","Upstream",
     &"Downstream","u(0,h/2)","u(L,h/2)","Inlet","Outlet","SUMIN"
     &,"SUMOUT","SUMLOWER"'
	close(3)
	
      ELSE

	Call readRestart        
*
      END IF 

	IF (NQLS.EQ.1) CALL SAMPI2

	endsubroutine
!========================================================================
	SUBROUTINE Readrestart

	Include 'common.txt'
	include 'property.txt'

	WRITE (*,*) ' READ THE RESTART FILE'
      OPEN (4,FILE='DSMC2.RES',STATUS='OLD',FORM='UNFORMATTED')
!      READ (4) ALPI,ALPN,ALPT,APX,APY,BME,BMEJ,BMR,BMRJ,BOLTZ,CB,
!     &           CC,CCG,CG,CH,COL,CS,CSR,CSS,CT,CW,CWRX,CWRY,DTM,FNDJ,
!     &           FNUM,FSPJ,FTMP,FVJ,FH,FW,IB,IC,IFCX,IFCY,IIS,IJET,IPL,
!     &           IPS,IR,ISC,ISCG,ISG,ISP,ISPR,ISURF,LFLX,LFLY,LIMJ,LIMS,
!     &           MOVT,NCOL,NCX,NCY,NIS,NM,NPS,NSCX,NSCY,NSMP,NPR,NPT,
!     &           NSP,PI,PP,PR,PV,RPX,RPY,SELT,SEPT,SP,SPI,SPM,SPR,TIME,
!     &           TIMI,TMPJ,TSURF,VFX,VFY,WJ,BMEINLET,BMEOUTLET,PIN,POUT,
!     &		   IRC(:,MNM),IPLOT	
!	READ(4) XGRID,YGRID,INZ,XYN

!from common.txt
	read(4)COL,MOVT,NCOL,SELT,SEPT,CS,CSR,CSS,NM,PP,PV,IPL,IPS,IR,
     &	IRC,PR,CC,CG,IC,ISC,CCG,ISCG,IG,NCX,NCY,IFCX,IFCY,CWRX,CWRY,
     &	APX,RPX,APY,RPY,SP,SPM,ISP,SPR,ISPR,CT,TIME,NPR,NSMP,
     &	FND,FTMP,TIMI,FSP,ISPD,VFX,VFY,
     &	FNUM,DTM,NIS,NSP,NPS,NPT,NSCX,NSCY,CB,IB,ISURF,LIMS,IIS,ISG,
     &	TSURF,IJET,LIMJ,TMPJ,FNDJ,FVJ,FSPJ,BMEJ,BMRJ,BME,BMR,CW,
     &	FW,CH,FH,WJ,LFLX,LFLY,ALPI,ALPN,ALPT,XGRID,YGRID,INZ,XYN,
     &	X_B,Y_B,PI,SPI,BOLTZ,
     &	BMEinlet,BMEoutlet,BMElower,IPROB,IPLOT,NBX,NBY,BF,
     &	H_BUFFER,BF,NBX,NBY,
!from property.txt
     &	SS,POUT,PIN

	close(4)

	NPT=600000
	
	end subroutine
!========================================================================
	SUBROUTINE writerestart
	
	Include 'common.txt'
	include 'property.txt'
	
!	WRITE (*,*) ' WRITING RESTART AND OUTPUT FILES',NPR,'  OF ',NPT
      OPEN (4,FILE='DSMC2.RES',FORM='UNFORMATTED')
!      WRITE (4) ALPI,ALPN,ALPT,APX,APY,BME,BMEJ,BMR,BMRJ,BOLTZ,CB,
!     &          CC,CCG,CG,CH,COL,CS,CSR,CSS,CT,CW,CWRX,CWRY,DTM,FNDJ,
!     &          FNUM,FSPJ,FTMP,FVJ,FH,FW,IB,IC,IFCX,IFCY,IIS,IJET,IPL,
!     &          IPS,IR,ISC,ISCG,ISG,ISP,ISPR,ISURF,LFLX,LFLY,LIMJ,LIMS,
!     &          MOVT,NCOL,NCX,NCY,NIS,NM,NPS,NSCX,NSCY,NSMP,NPR,NPT,NSP,
!     &          PI,PP,PR,PV,RPX,RPY,SELT,SEPT,SP,SPI,SPM,SPR,TIME,TIMI,
!     &          TMPJ,TSURF,VFX,VFY,WJ,BMEINLET,BMEOUTLET,PIN,POUT,
!     &		  IRC(:,MNM),IPLOT	
!	WRITE(4) XGRID,YGRID,INZ,XYN
!      CLOSE (4)
!from common.txt
	write(4)COL,MOVT,NCOL,SELT,SEPT,CS,CSR,CSS,NM,PP,PV,IPL,IPS,IR,
     &	IRC,PR,CC,CG,IC,ISC,CCG,ISCG,IG,NCX,NCY,IFCX,IFCY,CWRX,CWRY,
     &	APX,RPX,APY,RPY,SP,SPM,ISP,SPR,ISPR,CT,TIME,NPR,NSMP,
     &	FND,FTMP,TIMI,FSP,ISPD,VFX,VFY,
     &	FNUM,DTM,NIS,NSP,NPS,NPT,NSCX,NSCY,CB,IB,ISURF,LIMS,IIS,ISG,
     &	TSURF,IJET,LIMJ,TMPJ,FNDJ,FVJ,FSPJ,BMEJ,BMRJ,BME,BMR,CW,
     &	FW,CH,FH,WJ,LFLX,LFLY,ALPI,ALPN,ALPT,XGRID,YGRID,INZ,XYN,
     &	X_B,Y_B,PI,SPI,BOLTZ,
     &	BMEinlet,BMEoutlet,BMElower,IPROB,IPLOT,NBX,NBY,BF,
     &	H_BUFFER,BF,NBX,NBY,
!from property.txt
     &	SS,POUT,PIN

	close(4)



	end subroutine
!========================================================================
      SUBROUTINE MESH_SET
	
	Include 'common.txt'
      
	FW=CB(2)-CB(1)
      FH=CB(4)-CB(3)
      CG(1,1)=CB(1)
      IF (IFCX.EQ.0) THEN
        CW=FW/NCX
*--CW is the uniform cell width
      ELSE
        RPX=CWRX**(1./(NCX-1.))
*--RPX is the ratio in the geometric progression
        APX=(1.-RPX)/(1.-RPX**NCX)
*--AP is the first term of the progression
      END IF
      CG(4,1)=CB(3)
      IF (IFCY.EQ.0) THEN
        CH=FH/NCY
*--CH is the uniform cell height
      ELSE
        RPY=CWRY**(1./(NCY-1.))
*--RPY is the ratio in the geometric progression
        APY=(1.-RPY)/(1.-RPY**NCY)
*--APY is the first term of the progression
      END IF
      IF (IFCX.EQ.1) THEN
        APX=(1.-RPX)/APX
        RPX=LOG(RPX)
*--APX and RPX are now the convenient terms in eqn (12.1)
      END IF
      IF (IFCY.EQ.1) THEN
        APY=(1.-RPY)/APY
        RPY=LOG(RPY)
*--APY and RPY are now the convenient terms in eqn (12.1)
      END IF

	DO 500 MY=1,NCY
        DO 450 MX=1,NCX
          M=(MY-1)*NCX+MX
*--M is the cell number
          CT(M)=FTMP
*--the macroscopic temperature is set to the freestream temperature
*--set the x coordinates
          IF (MX.EQ.1) CG(1,M)=CG(1,1)
          IF (MX.GT.1) CG(1,M)=CG(2,M-1)
          IF (IFCX.EQ.0) THEN
            CG(2,M)=CG(1,M)+CW
          ELSE
            CG(2,M)=CG(1,M)+FW*APX*RPX**(MX-1)
          END IF
          CG(3,M)=CG(2,M)-CG(1,M)
*--set the y coordinates
          IF (MY.EQ.1) CG(4,M)=CG(4,1)
          IF (MY.GT.1.AND.MX.EQ.1) CG(4,M)=CG(5,M-1)
          IF (MY.GT.1.AND.MX.GT.1) CG(4,M)=CG(4,M-1)
          IF (IFCY.EQ.0) THEN
            CG(5,M)=CG(4,M)+CH
          ELSE
            CG(5,M)=CG(4,M)+FH*APY*RPY**(MY-1)
          END IF
          CG(6,M)=CG(5,M)-CG(4,M)
          CC(M)=CG(3,M)*CG(6,M)
          DO 420 L=1,MNSG
            DO 410 K=1,MNSG
              CCG(2,M,L,K)=RF(0)
              CCG(1,M,L,K)=SPM(1,1,1)*300.*SQRT(FTMP/300.)
410         CONTINUE
420       CONTINUE
*--the maximum value of the (rel. speed)*(cross-section) is set to a
*--reasonable, but low, initial value and will be increased as necessary
450     CONTINUE
500   CONTINUE

	END SUBROUTINE 
!========================================================================
	SUBROUTINE GRID_NOZZLE
	Include 'common.txt'

	INTEGER::I,J
	REAL::M,X,Y,DX,DY_CONV,DY_DIV,DY_throat,ythroat(ncy+1)


	DX=(cb(2)-cb(1))/NCX	
	DY_CONV=(cb(4)-cb(6))/ncy
	DY_DIV=(cb(4)-cb(3))/ncy

	DY_throat=(cb(4)-cb(5))/NCY
	
	xgrid(:,:)=0.0
	ygrid(:,:)=0.0

	ythroat(1)=cb(5)
	do j=2,ncy+1
		ythroat(j)=ythroat(j-1)+dy_throat
	enddo

	xgrid(1,1)=cb(1)
	do j=1,ncy+1
	do i=2,ncx+1
		xgrid(i,j)=xgrid(i-1,j)+dx
	enddo
	enddo

	ygrid(1:inz(1+1),1)=cb(6)
	DO J=1, NCY+1
		!the region before nozzle==>
			DO I=1, inz(1)+1
			  if(j>1) ygrid(i,j)=ygrid(i,j-1)+DY_CONV
			ENDDO
		
		!the converging region of the nozzle==>
		do i=inz(1)+2,inz(2)
*		 call LinearInterpolation(xgrid(inz(1)+1,j),ygrid(inz(1)+1,j),
*     &			xgrid(inz(2)+1,j),ythroat(j),xgrid(i,j),ygrid(i,j))
		  XP1=xgrid(inz(1)+1,j)
		  YP1=ygrid(inz(1)+1,j)
		  XP2=xgrid(inz(2)+1,j)
		  YP2=ythroat(j)
		  m=(YP2-YP1)/(XP2-XP1)
		  ygrid(i,j)=m*(xgrid(i,j)-XP1)+YP1

		enddo
		XYN(1,1)=Xgrid(inz(1)+1,1)
		XYN(1,2)=YGRID(inz(1)+1,1)
		ygrid(inz(2)+1,j)=ythroat(j)

		XYN(2,1)=XGRID(INZ(2)+1,1)
		XYN(2,2)=YGRID(INZ(2)+1,1)


		!the region after the nozzle==>
			DO I=inz(3)+1, ncx+1
				!if(j>1)	ygrid(i,j)=ygrid(i,j-1)+dy ! Diverging part
				if(j==1)ygrid(i,1)=CB(3)	!cb(5)	!FULLY CONVERGENT
				if(j>1)	ygrid(i,j)=ygrid(i,j-1)+(CB(4)-CB(3))/ncy !CB(5) FULLY CONVERGENT
			ENDDO
		XYN(3,1)=XGRID(INZ(3)+1,1)
		XYN(3,2)=YGRID(INZ(3)+1,1)

		!the diverging region of the nozzle==>
		do i=inz(2)+2,inz(3)
*		 call LinearInterpolation(xgrid(inz(2)+1,j),ygrid(inz(2)+1,j),
*     &		xgrid(inz(3)+1,j),ygrid(inz(3)+1,j),xgrid(i,j),ygrid(i,j))
			  XP1=xgrid(inz(2)+1,j)
			  YP1=ygrid(inz(2)+1,j)
			  XP2=xgrid(inz(3)+1,j)
			  YP2=ygrid(inz(3)+1,j)
			  m=(YP2-YP1)/(XP2-XP1)
			  ygrid(i,j)=m*(xgrid(i,j)-XP1)+YP1

		enddo
	ENDDO

!	BUFFER REGION
	CW=(CB(2)-CB(1))/NCX
	CH=(CB(4)-CB(3))/NCY

	Y_B(1,1)=CB(9)
	DO J=1,NBY+1
		X_B(1,J)=CB(7)
		IF(J.GT.1)THEN
			Y_B(1,J)=Y_B(1,J-1)+CH
		ENDIF
		DO I=2,NBX+1
			X_B(I,J)=X_B(I-1,J)+CW
			Y_B(I,J)=Y_B(1,J)
		ENDDO
	ENDDO

	ENDSUBROUTINE GRID_NOZZLE
!========================================================================
	SUBROUTINE PLOT_NOZZLE_2_ZONE
	
	Include 'common.txt'


	INTEGER::I,J

!	OPEN(1,FILE='GRID.PLT')
!	WRITE (1,*)'Variables=X,Y'
!	WRITE (1,*)'ZONE   I=',NCY,',  J=',NCX,',  F=POINT'
!	DO I=1,NCX
!	DO J=1,NCY
!		WRITE(1,*)XGRID(I,J),YGRID(I,J)
!	ENDDO
!	ENDDO
!	CLOSE(1)

	OPEN (7,FILE='GRID_2_ZONE.PLT')
	WRITE(7,*)'TITLE = "Subsonic"'
      WRITE (7,*)'Variables=X,Y'
	WRITE(7,*)'ZONE T="MAIN ZONE", I=',NCX+1,', J=',NCY+1,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	DO J=1,NCY+1
		DO I=1,NCX+1
			WRITE (7,'(2ES20.4)') XGRID(I,J),YGRID(I,J)
		ENDDO
	END DO

	WRITE(7,*)'ZONE T="BUFFER ZONE", I=',NBX+1,', J=',NBY+1,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	DO J=1,NBY+1
		DO I=1,NBX+1
			WRITE (7,'(2ES20.4)') X_B(I,J),Y_B(I,J)
		ENDDO
	END DO

      CLOSE (7)
	
	ENDSUBROUTINE PLOT_NOZZLE_2_ZONE
!========================================================================
	SUBROUTINE PLOT_GRID_2_ZONE
      include 'common.txt'
	include 'property.txt' 



	OPEN (4,FILE='GRID.PLT')
	WRITE(4,*)'TITLE = "Subsonic"'
      WRITE (4,*)'Variables=X,Y'
	WRITE(4,*)'ZONE T="MAIN ZONE", I=',NCX,', J=',NCY,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	DO N=1,MNC
	J=(N-0.1)/NCX+1
	I=N-(J-1)*NCX
	
!	IF((I.LE.LFLX).AND.(J.GE.LFLY))THEN
!	IF (IPROB(N).NE.0)THEN
	 CALL PROPERTIES (N)

       WRITE (4,'(2ES20.4)') XC,YC
!	ENDIF
!	ENDIF
	END DO

	WRITE(4,*)'ZONE T="BUFFER ZONE", I=',NBX,', J=',NBY,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	DO N=1,MNC
	J=(N-0.1)/NCX+1
	I=N-(J-1)*NCX
	
!	IF((I.GE.LFLX))THEN
!	IF (IPROB(N).NE.0)THEN
	 CALL PROPERTIES (N)

       WRITE (4,'(2ES20.4)') XC,YC
!	ENDIF
!	ENDIF
	END DO

      CLOSE (4)

	END SUBROUTINE
!========================================================================
	subroutine LinearInterpolation(x1,y1,x2,y2,x,y)
	implicit none
	real,intent(in)::x1,y1,x2,y2,x
	real,intent(out)::y
	real::m

	m=(y2-y1)/(x2-x1)
	y=m*(x-x1)+y1

	endsubroutine LinearInterpolation
!========================================================================
      SUBROUTINE MESH_SET_NOZZLE      
      
	Include 'common.txt'    
	


      CALL GRID_NOZZLE
	!CALL PLOT_NOZZLE_2_ZONE

	FW=CB(2)-CB(1)
      FH=CB(4)-CB(3)
      CG(1,1)=CB(1)
        CW=FW/NCX
*--CW is the uniform cell width
      CG(4,1)=CB(3)
        CH=FH/NCY
*--CH is the uniform cell height
      DO 500 MY=1,NCY
        DO 450 MX=1,NCX
          M=(MY-1)*NCX+MX
*--M is the cell number
          CT(M)=FTMP
*--the macroscopic temperature is set to the freestream temperature
*--set the x coordinates
		CG(1,M)=XGRID(MX,MY)	
		CG(2,M)=XGRID(MX+1,MY)
          CG(3,M)=CG(2,M)-CG(1,M)
*--set the y coordinates
          CG(4,M)=YGRID(MX,MY)		
		CG(5,M)=YGRID(MX+1,MY)
		CG(6,M)=YGRID(MX+1,MY+1)
		CG(7,M)=YGRID(MX,MY+1)
          CC(M)=CG(3,M)*(ABS(CG(7,M)-CG(4,M))+ABS(CG(6,M)-CG(5,M)))/2. !THE AREA OF THE CELL

          if((CG(6,M)<CG(5,M)).or.(CG(7,M)<CG(4,M))) print*,'ERROR  *5*'
		DO 420 L=1,MNSG
            DO 410 K=1,MNSG
              CCG(2,M,L,K)=RF(0)
              CCG(1,M,L,K)=SPM(1,1,1)*300.*SQRT(FTMP/300.)
410         CONTINUE
420       CONTINUE
*--the maximum value of the (rel. speed)*(cross-section) is set to a
*--reasonable, but low, initial value and will be increased as necessary
450     CONTINUE
500   CONTINUE

!	BUFFER REGION
	MNCM=NCX*NCY
	DO J=1,NBY
	  DO I=1,NBX
		M=MNCM+(J-1)*NBX+I
		CG(1,M)=X_B(I,J)
		CG(2,M)=X_B(I+1,J)
		CG(3,M)=CG(2,M)-CG(1,M)
		CG(4,M)=Y_B(I,J)
		CG(5,M)=Y_B(I+1,J)
		CG(6,M)=Y_B(I+1,J+1)
		CG(7,M)=Y_B(I,J+1)
          CC(M)=CG(3,M)*(ABS(CG(7,M)-CG(4,M))+ABS(CG(6,M)-CG(5,M)))/2. !THE AREA OF THE CELL
		DO 423 L=1,MNSG
            DO 413 K=1,MNSG
              CCG(2,M,L,K)=RF(0)
              CCG(1,M,L,K)=SPM(1,1,1)*300.*SQRT(FTMP/300.)
413         CONTINUE
423       CONTINUE
	  ENDDO
	ENDDO

*
*--set sub-cells
*
      DO 600 N=1,MNC
        DO 550 M=1,NSCY
          DO 520 K=1,NSCX
            L=(N-1)*NSCX*NSCY+(M-1)*NSCX+K
            ISC(L)=N
520       CONTINUE
550     CONTINUE
600   CONTINUE


	!determine each node's neighbors


	END SUBROUTINE 
!========================================================================
	SUBROUTINE ROWLocalSubCell(X1,Y1,X2,Y2,X3,Y3,X4,Y4,x,y,NROW)
	IMPLICIT NONE
	!this subroutine gets 4 points of a cell and then calculate the subcell 
	!number in which the molecule lies
	REAL,INTENT(IN)::X1,Y1,X2,Y2,X3,Y3,X4,Y4,X,y
	INTEGER,INTENT(OUT)::NROW
	REAL::YM1,YM2,YSC,M

	YM1=(Y1+Y4)/2.
	YM2=(Y2+Y3)/2.
*	CALL LinearInterpolation(X1,YM1,X2,YM2,X,YSC)

	m=(YM2-YM1)/(x2-x1)
	YSC=m*(x-x1)+YM1

	IF(Y<YSC) THEN
		NROW=1
	ELSE 
		NROW=2
	ENDIF

	ENDSUBROUTINE ROWLocalSubCell
!========================================================================
*   INIT2.FOR
	SUBROUTINE INIT2
*
*--initialization subroutine
*
      Include 'common.txt'

	REAL::MM,YLB,XP1,XP2,YP1,YP2,XN1,XN2,YN1,YN2,Y1,Y2

*
*--set constants
*
      PI=3.141592654
      SPI=SQRT(PI)
      BOLTZ=1.380622E-23
*
*--set data variables to default values that they retain if the data
*----does not reset them to specific values
      FND=0.
      FTMP=273.
      VFX=0.
      VFY=0.
      IFCX=0
      IFCY=0
      LFLX=0
      LFLY=0
      ALPI(1)=-1.
      ALPI(2)=-1.
	alpi(3)=-1.
      DO 100 N=1,4
        IB(N)=3
        DO 50 L=1,MNSP
          ISP(L)=1
          FSP(L)=0.
          BME(N,L)=0.
          BMR(N,L)=0.
50      CONTINUE
100   CONTINUE
      DO 200 L=1,MNSP
        BMEJ(L)=0.
        BMRJ(L)=0.
200   CONTINUE
*
      CALL DATA2
*
*--set additional data on the gas
*]
      IF (MNSP.EQ.1) ISPD=0
      DO 300 N=1,MNSP
        DO 250 M=1,MNSP
          IF ((ISPR(3,N).EQ.0).AND.(M.NE.N)) THEN
            SPR(1,N,M)=SPR(1,N,N)
            SPR(2,N,M)=SPR(2,N,N)
            SPR(3,N,M)=SPR(3,N,N)
          END IF
          IF ((ISPD.EQ.0).OR.(N.EQ.M)) THEN
            SPM(1,N,M)=0.25*PI*(SP(1,N)+SP(1,M))**2
*--the collision cross section is assumed to be given by eqn (1.35)
            SPM(2,N,M)=0.5*(SP(2,N)+SP(2,M))
            SPM(3,N,M)=0.5*(SP(3,N)+SP(3,M))
            SPM(4,N,M)=0.5*(SP(4,N)+SP(4,M))
*--mean values are used for ISPD=0
          ELSE
            SPM(1,N,M)=PI*SPM(1,N,M)**2
*--the cross-collision diameter is converted to the cross-section
          END IF
          SPM(5,N,M)=(SP(5,N)/(SP(5,N)+SP(5,M)))*SP(5,M)
*--the reduced mass is defined in eqn (2.7)
          SPM(6,N,M)=GAM(2.5-SPM(3,N,M))
250     CONTINUE
300   CONTINUE
*
*--initialise variables
*
      TIME=0.
      NM=0
      NPR=0
      NCOL=0
      MOVT=0.
      SELT=0.
      SEPT=0.
*
      DO 400 M=1,MNSP
        DO 350 N=1,MNSP
          COL(M,N)=0.
350     CONTINUE
400   CONTINUE

*	CALL MESH_SET
	CALL MESH_SET_NOZZLE
*
*
*--set sub-cells
*
      DO 600 N=1,MNC
        DO 550 M=1,NSCY
          DO 520 K=1,NSCX
            L=(N-1)*NSCX*NSCY+(M-1)*NSCX+K
            ISC(L)=N
520       CONTINUE
550     CONTINUE
600   CONTINUE
*
      IF (IIS.GT.0.AND.ISG.GT.0) THEN
*--if IIS=1 generate initial gas with temperature FTMP
*
        DO 650 L=1,MNSP
          REM=0
          IF (IIS.EQ.1) VMP=SQRT(2.*BOLTZ*FTMP/SP(5,L))
*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
          DO 620 N=1,MNC
            IPROB(N)=1
            IF (LFLX.NE.0) THEN
              NY=(N-1)/NCX+1
              NX=N-(NY-1)*NCX
* NX and NY are the cell column and row
!              IF ((LFLX.GT.0.AND.LFLY.GT.0).AND.
!     &            (NX.LT.LFLX.AND.NY.LT.LFLY)) IPROB(N)=0
!              IF ((LFLX.GT.0.AND.LFLY.LT.0).AND.
!     &            (NX.LT.LFLX.AND.NY.GT.-LFLY)) IPROB(N)=0
!              IF ((LFLX.LT.0.AND.LFLY.GT.0).AND.
!     &            (NX.GT.-LFLX.AND.NY.LT.LFLY)) IPROB(N)=0
!              IF ((LFLX.LT.0.AND.LFLY.LT.0).AND.
!     &            (NX.GT.-LFLX.AND.NY.GT.-LFLY)) IPROB(N)=0
            END IF
            IF (IPROB(N).EQ.1) THEN
              A=FND*CC(N)*FSP(L)/FNUM+REM
*--A is the number of simulated molecules of species L in cell N to
*--simulate the required concentrations at a total number density of FND
              IF (N.LT.MNC) THEN
                MM=A
                REM=(A-MM)
*--the remainder REM is carried forward to the next cell
              ELSE
                MM=NINT(A)
              END IF
              IF (MM.GT.0) THEN
                DO 604 M=1,MM
                  IF (NM.LT.MNM) THEN
*--round-off error could have taken NM to MNM+1
                    NM=NM+1
                    IPS(NM)=L
                    PP(1,NM)=CG(1,N)+RF(0)*(CG(2,N)-CG(1,N))
                    NCOLM=(PP(1,NM)-CG(1,N))*(NSCX-.001)/CG(3,N)+1
!                   PP(2,NM)=CG(4,N)+RF(0)*(CG(5,N)-CG(4,N))

*     			  CALL LinearInterpolation(CG(1,N),CG(4,N),CG(2,N) 
*     &			 ,CG(5,N),pp(1,nm),y1)
				  XP1=CG(1,N)
				  YP1=CG(4,N)
				  XP2=CG(2,N)
				  YP2=CG(5,N)
				  MM=(YP2-YP1)/(XP2-XP1)
				  Y1=MM*(PP(1,NM)-XP1)+YP1
				  
*				  CALL LinearInterpolation(CG(1,N),CG(7,N),CG(2,N) 
*     &			 ,CG(6,N),pp(1,nm),y2)
				  XN1=CG(1,N)
				  YN1=CG(7,N)
				  XN2=CG(2,N)
				  YN2=CG(6,N)
				  MM=(YN2-YN1)/(XN2-XN1)
				  Y2=MM*(PP(1,NM)-XN1)+YN1

				  CALL Y_LOWER_BOUND(PP(1,NM),YLB)
				  IF(nm==484)THEN
					YLB=YLB+0.0
				  ENDIF
				  if(y2<y1) print*,'ERROR  *3*'
				  PP(2,NM)=y1+RF(0)*(y2-y1)
				  IF(PP(2,NM)<YLB)THEN
					PRINT*,'**ERROR INIT ** Y < YLB'
				  CALL Y_LOWER_BOUND(PP(1,NM),YLB)
				  ENDIF

				  !NROW=(PP(2,NM)-CG(4,N))*(NSCY-.001)/CG(6,N)+1

!	              IRC(2,NM)=(N-IRC(1,nm))/ncx+1  !(PP(2,NM)+1E-15)/CG(6,N) + 1
                    !NROW=(PP(2,NM)-CG(4,N))*(NSCY-.001)/CG(6,N)+1
		  		  CALL ROWLocalSubCell(cg(1,N),cg(4,N),cg(2,N),cg(5,N),
     &					cg(2,N),cg(6,N),cg(1,N),cg(7,N),
     &					PP(1,NM),pp(2,nm),NROW)


                    IPL(NM)=(N-1)*NSCX*NSCY+(NROW-1)*NSCX+NCOLM
*--species, position, and sub-cell number have been set
                    DO 602 K=1,3
                    CALL RVELC(PV(K,NM),A,VMP)
602                 CONTINUE
                    PV(1,NM)=PV(1,NM)+VFX
                    PV(2,NM)=PV(2,NM)+VFY
*--velocity components have been set
*--set the rotational energy
                    IF (ISPR(1,L).GT.0) CALL SROT(PR(NM),FTMP,ISPR(1,L))
                  END IF
604             CONTINUE
              END IF
            END IF
620       CONTINUE
650     CONTINUE

      IC(2,:,:)=0
	DO  N=1,NM
        MSC=IPL(N)
        MC=ISC(MSC)
        IC(2,MC,1)=IC(2,MC,1)+1
	  IF(IC(2,MC,1)==0) PRINT*,'Error M **********'
	ENDDO

        WRITE (*,99001) NM
99001   FORMAT (' ',I6,' MOLECULES')
      END IF      
*--now calculate the number that enter in jet
      IF (IJET.GT.0) THEN
*--molecules enter from a jet
        DO 750 L=1,MNSP
          VMP=SQRT(2.*BOLTZ*TMPJ/SP(5,L))
*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
          SC=FVJ/VMP
*--SC is the inward directed speed ratio
          IF (ABS(SC).LT.10.1) A=(EXP(-SC*SC)+SPI*SC*(1.+ERF(SC)))
     &                           /(2.*SPI)
          IF (SC.GT.10.) A=SC
          IF (SC.LT.-10.) A=0.
*--A is the non-dimensional flux of eqn (4.22)
          J2=LIMJ(2)
          J3=LIMJ(3)
          IF (IJET.EQ.1.OR.IJET.EQ.2) THEN
            WJ=CG(2,J3)-CG(1,J2)
*--WJ is the width of the jet
          ELSE
            WJ=CG(5,(J3-1)*NCX+1)-CG(4,(J2-1)*NCX+1)
          END IF
          BMEJ(L)=FNDJ*FSPJ(L)*A*VMP*DTM*WJ/FNUM
!          WRITE (*,*) ' entering mols in jet ',BMEJ(L)
*--the entering number is given by eqn (12.5)
750     CONTINUE
      END IF
      RETURN
      END
!======================================================================
	SUBROUTINE Y_LOWER_BOUND(X,YLB)

      Include 'common.txt'
	
	REAL, INTENT(IN)::X
	REAL, INTENT(OUT)::YLB
	REAL M,Y1,Y2,X1,X2

	IF(X.LE.XYN(1,1))THEN
	YLB=CB(6)	!6 BECAUSE THE INLET ZONE IS HIGHER  !CB(3)
	ELSEIF((X.GT.XYN(1,1)).AND.(X.LT.XYN(2,1)))THEN
*	CALL LinearInterpolation(XYN(1,1),XYN(1,2),
*    &	XYN(2,1),XYN(2,2),X,YLB)
	X1=XYN(1,1)
	X2=XYN(2,1)
	Y1=XYN(1,2)
	Y2=XYN(2,2)
	m=(y2-y1)/(x2-x1)
	YLB=m*(x-x1)+y1

	ELSEIF((X.GT.XYN(2,1)).AND.(X.LT.XYN(3,1)))THEN
*	ELSEIF((X.GT.XYN(2,1)))THEN	!fully convergent
*	YLB=CB(5)	!fully convergent
	X1=XYN(2,1)
	X2=XYN(3,1)
	Y1=XYN(2,2)
	Y2=XYN(3,2)
	m=(y2-y1)/(x2-x1)
	YLB=m*(x-x1)+y1

	ELSEIF((X.GT.XYN(3,1)).AND.(X.LT.CB(2)))THEN
	YLB=CB(3)
	ELSEIF(X.GT.CB(2))THEN
	YLB=CB(9)
	ENDIF

	END SUBROUTINE Y_LOWER_BOUND
!======================================================================
      SUBROUTINE MOVE2
*
*--the NM molecules are moved over the time interval DTM
*
      Include 'common.txt'
*
      integer::nc,nx,ncnew,nxnew,mec=0
	REAL X,XI,Y,YI,YLB,DX,DY,YS,X1,X2,Y1,Y2,XD	!,XC	!
	REAL XSU,XSD,M,M1,M2,alfa,VT,VN
	IFT=-1
*--a negative IFT indicates that molecules have not entered at this step
      N=0
	CW=(CB(2)-CB(1))/NCX
	CH=(CB(4)-CB(3))/NCY

100   N=N+1
      IF (N.LE.NM) THEN
        IF (IFT.LT.0) AT=DTM
        IF (IFT.GT.0) AT=RF(0)*DTM
*--the time step is a random fraction of DTM for entering molecules
150     MOVT=MOVT+1
        MSC=IPL(N)
        MC=ISC(MSC)
	  X=PP(1,N)
	  Y=PP(2,N)
        IF(X.LT.CB(2))THEN
		IF (IFCX.EQ.0) THEN
            MCX=FLOOR((X-CB(1))/CW)+1
		  IF(X<CG(1,MCX)) MCX=MCX-1
*		  IF((floor(x/cb(2)*real(ncx))+1).ne.mcx) PRINT*, '**ERROR 5**'
          ELSE
            XD=(X-CB(1))/FW+1.E-6
            MCX=1.+(LOG(1.-XD*APX))/RPX
*--the cell number is calculated from eqn (12.1)
          END IF
          IF (MCX.LT.1) MCX=1
          IF (MCX.GT.NCX) MCX=NCX
*--MCX is the new cell column (note avoidance of round-off error)
          IF (IFCY.EQ.0) THEN
		  CALL Y_LOWER_BOUND(X,YLB)
		  MCY=(Y-YLB)/(CB(4)-YLB)*NCY+0.99999	!(Y-CB(3))/CH+0.99999
*		  MCY=(Y-CB(3))/CH+0.99999
          ELSE
            YD=(Y-CB(3))/FH+1.E-6
            MCY=1.+(LOG(1.-YD*APY))/RPY
*--the cell number is calculated from eqn (12.1)
          END IF
          IF (MCY.LT.1) MCY=1
          IF (MCY.GT.NCY) MCY=NCY
*--MCY is the new cell row (note avoidance of round-off error)
          MC2=(MCY-1)*NCX+MCX
		
	  elseif(x.gt.cb(2))then
		mcx=floor((x-cb(2))/cw)+1
	    if(x.lt.cg(1,ncx*ncy+mcx)) mcx=mcx-1
		if(mcx.lt.1) mcx=1
	    if(x.gt.cg(2,ncx*ncy+mcx)) mcx=mcx+1
		if(mcx.gt.nbx) mcx=nbx
		
		mcy=floor((y-cb(9))/ch)+1
	    if(y.lt.cg(4,ncx*ncy+mcx)) mcy=mcy-1
		if(mcy.lt.1) mcy=1
		if(y.gt.cg(6,ncx*ncy+(mcy-1)*nbx+mcx)) mcy=mcy+1
		if(mcy.gt.nby) mcy=nby
	    mc2=ncx*ncy+(mcy-1)*nbx+mcx
	  endif

	IF(MC2.NE.MC)THEN
	!PRINT*,'**ERROR**CELL NUMBER MOVE**',MC,MC2
	ENDIF
	!MC=MC2
*--MC is the initial cell number
        XI=PP(1,N)
        IF ((XI+0.00001*CG(3,1)).LT.CB(1).OR.
     &   (XI-0.00001*CG(3,MNC)).GT.CB(8)) THEN
          WRITE (*,*) ' MOL ',N,' X COORD OUTSIDE FLOW ',XI
          CALL REMOVE(N)
          GO TO 100
        END IF
        YI=PP(2,N)
*        IF ((YI+0.00001*(CG(6,1)-CG(5,1))).LT.CB(3).OR.
*     &   (YI-0.00001*(CG(6,1)-CG(5,1))).GT.CB(4)) THEN
        CALL Y_LOWER_BOUND(XI,YLB)
	  IF ((YI+0.00001*(CG(6,1)-CG(5,1))).LT.YLB.OR.
     &   (YI-0.00001*(CG(6,1)-CG(5,1))).GT.CB(4)) THEN
          WRITE (*,*) ' MOL ',N,' Y COORD OUTSIDE FLOW ',xi,YI,ylb
          CALL REMOVE(N)
          GO TO 100
        END IF
        DX=PV(1,N)*AT
        DY=PV(2,N)*AT
        X=XI+DX
        Y=YI+DY

        DO 200 KS=1,3
*--check the surfaces
          IF (ISURF(KS).GT.0) THEN
            IF (ISURF(KS).EQ.1.OR.ISURF(KS).EQ.2) THEN
              L1=LIMS(KS,1)
              IF (L1.LE.NCY) THEN
			  CALL Y_LOWER_BOUND(XI,YLB)
			  CALL Y_LOWER_BOUND(X,YS)	!YS=CG(4,(L1-1)*NCX+1)
              ELSE
                YS=CB(4)
                L1=L1-1
              END IF
              IF ((ISURF(KS).EQ.1.AND.(YI.GT.YLB.AND.Y.LT.YS)).OR.
     &            (ISURF(KS).EQ.2.AND.(YI.LT.YS.AND.Y.GT.YS))) THEN
                IF(ISURF(KS).EQ.1)THEN
				  IF(XI<XYN(1,1))THEN
					  XC=XI+(YS-YI)*DX/DY
					  IF(XC.GT.XYN(1,1).AND.X.GT.XYN(1,1))THEN
						X1=XYN(1,1);Y1=XYN(1,2)
						X2=XYN(2,1);Y2=XYN(2,2)
						M1=(Y2-Y1)/(X2-X1)
						M2=(Y-YI)/(X-XI)
						XC=(M1*X1-M2*XI+YI-Y1)/(M1-M2)
					  ENDIF
				  ELSEIF((XI>XYN(1,1)).AND.(XI<XYN(2,1)))THEN
					X1=XYN(1,1);Y1=XYN(1,2)
					X2=XYN(2,1);Y2=XYN(2,2)
					M1=(Y2-Y1)/(X2-X1)
					M2=(Y-YI)/(X-XI)
					XC=(M1*X1-M2*XI+YI-Y1)/(M1-M2)
					IF(XC.GT.XYN(2,1).AND.X.GT.XYN(2,1))THEN
						X1=XYN(2,1);Y1=XYN(2,2)
						X2=XYN(3,1);Y2=XYN(3,2)
						M1=(Y2-Y1)/(X2-X1)
						M2=(Y-YI)/(X-XI)
						XC=(M1*X1-M2*XI+YI-Y1)/(M1-M2)
					ENDIF
					IF(((XI<X).AND.((XC<XI).OR.(XC>X)))
     &					.OR.((XI>X).AND.((XC>XI).OR.(XC<X))))THEN
						PRINT*, 'ERROR IN XC'
					ENDIF
					!XC=(YI-Y1+(Y2-Y1)/(X2-X1)*X1-(Y-YI)/(X-XI)*XI)/
!     &				!	((Y2-Y1)/(X2-X1)-(Y-YI)/(X-XI))
*				  ELSEIF((XI>XYN(2,1)).AND.(XI<XYN(3,1)))THEN
				  ELSEIF(XI.GT.XYN(2,1).AND.XI.LT.XYN(3,1))THEN	
*					  XC=XI+(YS-YI)*DX/DY	!fully convergent
					X1=XYN(2,1);Y1=XYN(2,2)
					X2=XYN(3,1);Y2=XYN(3,2)
					M1=(Y2-Y1)/(X2-X1)
					M2=(Y-YI)/(X-XI)
					XC=(M1*X1-M2*XI+YI-Y1)/(M1-M2)
					IF((XC.GT.XYN(3,1).OR.(XC.LT.XI.AND.DX.GT.0))
     &						.AND.(XYN(3,1).LT.CB(2)))THEN
					    XC=XI+(YS-YI)*DX/DY
					ENDIF
					IF(((XI<X).AND.((XC<XI).OR.(XC>X)))
     &					.OR.((XI>X).AND.((XC>XI).OR.(XC<X))))THEN
						PRINT*, 'ERROR IN XC ==> DIVERGENT SECTION',XC
					ENDIF
				  ELSEIF((XI>XYN(3,1)).AND.(XI<CB(2)))THEN
					  XC=XI+(YS-YI)*DX/DY
				  ELSEIF(XI.GT.CB(2))THEN
					  GOTO 153
				  ENDIF

	          ELSEIF(ISURF(KS).EQ.2)THEN
				  XC=XI+(YS-YI)*DX/DY
			  ENDIF
                IF (XC.LE.CB(1).AND.IB(1).EQ.2) THEN !IB(1).EQ.2 SYMMETRY BOUNDARY
                  XC=2.*CB(1)-XC
                  PV(1,N)=-PV(1,N)
                END IF
                IF (XC.GE.CB(8).AND.IB(2).EQ.2) THEN !IB(2).EQ.2 SYMMETRY BOUNDARY
                  XC=2.*CB(8)-XC
                  PV(1,N)=-PV(1,N)
                END IF
                L2=LIMS(KS,2)
                L3=LIMS(KS,3)
                XSU=CG(1,L2)
                XSD=CG(2,L3)
                IF (XC.GT.XSU.AND.XC.LT.XSD) THEN  !WALL BOUNDARY
*--molecule collides with surface at XC
                  IF (IFCX.EQ.0) THEN
                    MC=FLOOR((XC-CB(1))/CW)+1
				  IF(XC<CG(1,MC)) MC=MC-1
                  ELSE
                    XD=(XC-CB(1))/FW+1.E-6
                    MC=1.+(LOG(1.-XD*APX))/RPX
*--the cell number is calculated from eqn (12.1)
                  END IF
                  IF (MC.LT.1) MC=1
                  IF (MC.GT.NCX) MC=NCX
                  MCS=MC-(L2-1)
                  IF (ISURF(KS).EQ.1) MC=MC+(L1-1)*NCX
                  IF (ISURF(KS).EQ.2)then
				   MC=MC+(L1-1)*NCX			!modified by Mirjalili
				   MCS=MC-L2
				ENDIF
*--MC is the cell number for the reflected molecule
                  IF (KS.EQ.2) MCS=MCS+LIMS(1,3)-LIMS(1,2)+1
*--MCS is the code number of the surface element
                  CALL Y_LOWER_BOUND(XC,YSC)
				AT=AT*(Y-YSC)/DY
				IF(MC.GT.NCX) PRINT*,'****ERROR*MC > NCX***'

                  CALL REFLECT2(N,KS,MCS,XC,YSC,MC)
                  GO TO 150
                END IF
              END IF
            END IF
153         IF (ISURF(KS).EQ.3.OR.ISURF(KS).EQ.4) THEN
              !L1=LIMS(KS,1)
              !IF (L1.LE.NCX) THEN
              !  XS=CG(1,L1)
              !ELSE
                XS=CB(2)
              !  L1=L1-1
              !END IF
              IF ((ISURF(KS).EQ.3.AND.(XI.GT.XS.AND.X.LT.XS)).OR.
     &            (ISURF(KS).EQ.4.AND.(XI.LT.XS.AND.X.GT.XS))) THEN
                YC=YI+(XS-XI)*DY/DX
                !L2=LIMS(KS,2)
                !L3=LIMS(KS,3)
                YSU=cb(9)	!CG(4,(L2-1)*NCX+1)
                YSD=cb(3)	!CG(5,(L3-1)*NCX+1)
                IF (YC.GT.YSU.AND.YC.LT.YSD) THEN
*--molecule collides with surface at YC
                  IF (IFCY.EQ.0) THEN
                    MC=(YC-cb(9))/CH+0.99999
                  ELSE
                    YD=(YC-CB(3))/FH+1.E-6
                    MC=1.+(LOG(1.-YD*APY))/RPY
*--the cell number is calculated from eqn (12.1)
                  END IF
                  IF (MC.LT.1) MC=1
                  IF (MC.GT.NCY) MC=nby-NCY
                  MCS=MC	!-(L2-1)
                  IF (ISURF(KS).EQ.3) MC=ncx*ncy+(mc-1)*nbx+1	!(MC-1)*NCX+L1
                  IF (ISURF(KS).EQ.4) MC=(MC-1)*NCX+L1-1
*--MC is the cell number for the reflected molecule
                  IF (KS.EQ.2) MCS=MCS+LIMS(1,3)-LIMS(1,2)+1
*--MCS is the code number of the surface element
                  AT=AT*(XS-XI)/DX
                  CALL REFLECT2(N,KS,MCS,XS,YC,MC)
                  GO TO 150
                END IF
              END IF
            END IF
          END IF
200     CONTINUE
        IF (X.LT.CB(1).OR.X.GT.cb(8)) THEN
          IF (X.LT.CB(1)) K=1
          IF (X.GT.cb(8)) K=2
*--intersection with boundary K
          IF (IB(K).EQ.2) THEN
*--specular reflection from the boundary (eqn (11.7))
            X=2.*CB(K)-X
            PV(1,N)=-PV(1,N)
          ELSE
*--molecule leaves flow
            CALL REMOVE(N)
            GO TO 100
          END IF
        END IF
        CALL Y_LOWER_BOUND(X,YS)
	  IF (Y.LT.YS.OR.Y.GT.CB(4)) THEN
          IF ((Y.LT.YS).and.(x.lt.cb(2))) K=3
		if ((Y.LT.YS).and.(x.gt.cb(2))) K=6
          IF (Y.GT.CB(4)) K=4
*--intersection with boundary K
          IF (IB(K).EQ.2) THEN
*--specular reflection from the boundary (eqn (11.7))
		  IF(K==3)THEN 
		      Y=2.*YS-Y
		      IF(X<XYN(1,1))THEN
				PV(2,N)=-PV(2,N)
			  ELSEIF(X.GT.XYN(1,1).AND.X.LT.XYN(2,1))THEN
*			  ELSEIF(X>XYN(1,1))THEN	!FULLY CONVERGENT 
				X1=XYN(1,1);Y1=XYN(1,2)
      			X2=XYN(2,1);Y2=XYN(2,2)
				M1=(Y2-Y1)/(X2-X1)
				ALFA=ATAN(M1)
				VT=PV(1,N)*COS(ALFA)+PV(2,N)*SIN(ALFA)
				VN=-PV(1,N)*SIN(ALFA)+PV(2,N)*COS(ALFA)
				PV(1,N)=VT*COS(ALFA)+VN*SIN(ALFA)
				PV(2,N)=VT*SIN(ALFA)-VN*COS(ALFA)
			  ELSEIF(X.GT.XYN(2,1).AND.X.LT.XYN(3,1))THEN
				X1=XYN(2,1);Y1=XYN(2,2)
      			X2=XYN(3,1);Y2=XYN(3,2)
				M1=(Y2-Y1)/(X2-X1)
				ALFA=ATAN(M1)
				VT=PV(1,N)*COS(ALFA)+PV(2,N)*SIN(ALFA)
				VN=-PV(1,N)*SIN(ALFA)+PV(2,N)*COS(ALFA)
				PV(1,N)=VT*COS(ALFA)+VN*SIN(ALFA)
				PV(2,N)=VT*SIN(ALFA)-VN*COS(ALFA)
			  ELSEIF(X.GT.XYN(3,1))THEN
				PV(2,N)=-PV(2,N)
			  ENDIF
	      ELSEIF(K==4)THEN
			  Y=2.*CB(K)-Y
			  PV(2,N)=-PV(2,N)
		  ENDIF
          ELSE
*--molecule leaves flow
            CALL REMOVE(N)
            GO TO 100
          END IF
        END IF
*
*      CALL LinearInterpolation(CG(1,MC),CG(4,MC),CG(2,MC),CG(5,MC),X,YL)
	  XP1=CG(1,MC)
	  YP1=CG(4,MC)
	  XP2=CG(2,MC)
	  YP2=CG(5,MC)
	  m=(YP2-YP1)/(XP2-XP1)
	  YL=m*(X-XP1)+YP1

*      CALL LinearInterpolation(CG(1,MC),CG(7,MC),CG(2,MC),CG(6,MC),X,YU)
	  XP1=CG(1,MC)
	  YP1=CG(7,MC)
	  XP2=CG(2,MC)
	  YP2=CG(6,MC)
	  m=(YP2-YP1)/(XP2-XP1)
	  YU=m*(X-XP1)+YP1

*        IF (X.LT.CG(1,MC).OR.X.GT.CG(2,MC).OR.Y.LT.CG(4,MC).OR.
*     &      Y.GT.CG(5,MC)) THEN
	  IF (X.LT.CG(1,MC).OR.X.GT.CG(2,MC).OR.Y.LT.YL.OR.Y.GT.YU) THEN
*--the molecule has moved from the initial cell
          if(x.lt.cb(2))then
		IF (IFCX.EQ.0) THEN
            MCX=FLOOR((X-CB(1))/CW)+1
		  IF(X<CG(1,MCX)) MCX=MCX-1

*		  IF((floor(x/cb(2)*real(ncx))+1).ne.mcx) PRINT*, '**ERROR 5**'
          ELSE
            XD=(X-CB(1))/FW+1.E-6
            MCX=1.+(LOG(1.-XD*APX))/RPX
*--the cell number is calculated from eqn (12.1)
          END IF
          IF (MCX.LT.1) MCX=1
          IF (MCX.GT.NCX) MCX=NCX
*--MCX is the new cell column (note avoidance of round-off error)
          IF (IFCY.EQ.0) THEN
		  CALL Y_LOWER_BOUND(X,YLB)
		  MCY=(Y-YLB)/(CB(4)-YLB)*NCY+0.99999	!(Y-CB(3))/CH+0.99999
*		  MCY=(Y-CB(3))/CH+0.99999
          ELSE
            YD=(Y-CB(3))/FH+1.E-6
            MCY=1.+(LOG(1.-YD*APY))/RPY
*--the cell number is calculated from eqn (12.1)
          END IF
          IF (MCY.LT.1) MCY=1
          IF (MCY.GT.NCY) MCY=NCY
*--MCY is the new cell row (note avoidance of round-off error)
          MC=(MCY-1)*NCX+MCX
	  elseif(x.gt.cb(2))then
		mcx=floor((x-cb(2))/cw)+1
	    if(x.lt.cg(1,ncx*ncy+mcx)) mcx=mcx-1
		if(mcx.lt.1) mcx=1
	    if(x.gt.cg(2,ncx*ncy+mcx)) mcx=mcx+1
		if(mcx.gt.nbx) mcx=nbx
		
		mcy=floor((y-cb(9))/ch)+1
	    if(y.lt.cg(4,ncx*ncy+mcx)) mcy=mcy-1
		if(mcy.lt.1) mcy=1
		if(y.gt.cg(6,ncx*ncy+(mcy-1)*nbx+mcx)) mcy=mcy+1
		if(mcy.gt.nby) mcy=nby
	    mc=ncx*ncy+(mcy-1)*nbx+mcx
	  endif
        END IF
        MSCX=((X-CG(1,MC))/CG(3,MC))*(NSCX-.001)+1
*       MSCY=((Y-CG(4,MC))/CG(6,MC))*(NSCY-.001)+1
        CALL ROWLocalSubCell(cg(1,MC),cg(4,MC),cg(2,MC),cg(5,MC),
     &	cg(2,MC),cg(6,MC),cg(1,MC),cg(7,MC),x,y,mscy)
	  MSC=(MSCY-1)*NSCX+MSCX+NSCX*NSCY*(MC-1)
*--MSC is the new sub-cell number
        IF (MSC.LT.1) MSC=1
        IF (MSC.GT.MNSC) MSC=MNSC
        IPL(N)=MSC
        PP(1,N)=X
	  CALL Y_LOWER_BOUND(X,YLB)
	  IF(Y<YLB)THEN
	    PRINT*,'*****ERROR** Y < YLB ***'
	  ENDIF
        PP(2,N)=Y
        GO TO 100
      ELSE IF (IFT.LT.0) THEN
        IFT=1
*--new molecules enter
        IF (NSMP.GE.1) CALL ENTER2
        N=N-1
        GO TO 100
      END IF
*	print*,"the number of exit molecules from a cell: ", mec
      RETURN
      END
!======================================================================
      SUBROUTINE ENTER2

*--new molecules enter at boundaries
	INCLUDE 'COMMON.TXT'
	INCLUDE 'PROPERTY.TXT'

	REAL::XP1,XP2,YP1,YP2,YL,YU,MSLOPE,YLB	

	IF (IIS.GT.0) THEN
*
*--calculate the number of molecules that enter at each time step
*--across the four sides of the simulated region
        DO 700 N=1,6
          IF (IB(N).EQ.1) THEN
!            WRITE (*,*) 'side',N
*--molecules enter from an external stream
            DO 660 L=1,MNSP
              VMP=SQRT(2.*BOLTZ*FTMP/SP(5,L))
*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
!              IF (N.EQ.1) SC=VFX/VMP
!              IF (N.EQ.2) SC=-VFX/VMP
              IF (N.EQ.3) SC=VFY/VMP
              IF (N.EQ.4) SC=-VFY/VMP
!              if (n.eq.5) sc=vfx/vmp
!              if (n.eq.6) SC=vfy/vmp
*--SC is the inward directed speed ratio
              IF (ABS(SC).LT.10.1) A=(EXP(-SC*SC)+SPI*SC*(1.+ERF(SC)))
     &                               /(2.*SPI)
              IF (SC.GT.10.) A=SC
              IF (SC.LT.-10.) A=0.
*--A is the non-dimensional flux of eqn (4.22)
!              IF (N.EQ.1.OR.N.EQ.2) THEN
!                BME(N,L)=FND*FSP(L)*A*VMP*DTM*FH/FNUM
!              ELSE
!                BME(N,L)=FND*FSP(L)*A*VMP*DTM*FW/FNUM
			if (n.GT.2) then
                BME(N,L)=FND*FSP(L)*A*VMP*DTM*FW/FNUM
              END IF              
!			WRITE (*,*) ' entering mols ',BME(N,L)
660         CONTINUE
          END IF
700     CONTINUE
      END IF

	do n=1,3
	L=1

	if(n.eq.1)then
	nlim=ncy
	elseif(n.eq.2)then
	nlim=nby   !out	
	elseif(n.eq.3)then
	nlim=nbx   !lower side
	endif

	Do NN=1,nlim
	if(n.eq.1) mc=(nn-1)*ncx+1
	if(n.eq.2) mc=ncx*ncy+nn*nbx  !out		 
	if(n.eq.3) mc=ncx*ncy+nn      !lower side	  

     	CALL PROPERTIES(mc)

	VMP=SQRT(2.*BOLTZ/SP(5,L)*temp)
	IF (N.EQ.1) SC=VEL(1)/VMP
	IF (N.EQ.2) SC=-VEL(1)/VMP
	if (n.eq.3) sc=-vel(2)/vmp

	IF (ABS(SC).LT.10.1) A=(EXP(-SC*SC)+SPI*SC*(1.+ERF(SC)))/(2.*SPI)
      IF (SC.GT.10.) A=SC
      IF (SC.LT.-10.) A=0.
	! for uniform cell ?
	
!	MM=(NN-1)*NCX+1	
	mm=mc
	if (N.eq.1) BMEinlet(NN)=DENN*A*VMP*DTM*(CG(7,MM)-CG(4,MM))/FNUM
	if (N.eq.2) BMEoutlet(NN)=DENN*A*VMP*DTM*(CG(6,MM)-CG(5,MM))/FNUM 
	if (n.eq.3) BMElower(nn)=denn*1*vmp*dtm*cg(3,mm)/fnum

	end do
	end do
*
      DO 100 N=1,6	!4
*--consider each boundary in turn
        IF (IB(N).EQ.1) THEN
          !IF (N.LT.3) NCS=NCY
          !IF (N.GT.2) NCS=NCX
		if(n.eq.1) ncs=ncy
		if(n.eq.2) ncs=nby
		if(n.eq.6) ncs=nbx

          DO 20 NC=1,NCS
		if (n.eq.1)	CALL PROPERTIES((nc-1)*ncx+1)	!.and.IPROB((NC-1)*ncx+1).NE.0)
		if (n.eq.2) CALL PROPERTIES(ncx*ncy+nc*nbx)	!.and.IPROB(NC*ncx).NE.0
		if (n.eq.6) call properties(ncx*ncy+nc)

            IF (LFLX.NE.0) THEN
*--bypass entry into the excluded region of the flow
              IF (N.EQ.1) THEN
                IF (LFLY.GT.0.AND.LFLX.GT.0.AND.NC.LT.LFLY) GO TO 20
                IF (LFLY.LT.0.AND.LFLX.GT.0.AND.NC.GT.LFLY) GO TO 20
              END IF
              IF (N.EQ.2) THEN
                IF (LFLY.GT.0.AND.LFLX.LT.0.AND.NC.LT.LFLY) GO TO 20
                IF (LFLY.LT.0.AND.LFLX.LT.0.AND.NC.GT.LFLY) GO TO 20
              END IF
              IF (N.EQ.3) THEN
                IF (LFLX.GT.0.AND.LFLY.GT.0.AND.NC.LT.LFLX) GO TO 20
                IF (LFLX.LT.0.AND.LFLY.GT.0.AND.NC.GT.LFLX) GO TO 20
              END IF
              IF (N.EQ.4) THEN
                IF (LFLX.GT.0.AND.LFLY.LT.0.AND.NC.LT.LFLX) GO TO 20
                IF (LFLX.LT.0.AND.LFLY.LT.0.AND.NC.GT.LFLX) GO TO 20
              END IF
            END IF

            DO 10 L=1,MNSP
*--consider each species in turn
              VMP=SQRT(2.*BOLTZ/SP(5,L)*FTMP)
			if ((N.LT.3).and.(n.eq.6)) VMP=SQRT(2.*BOLTZ/SP(5,L)*TEMP) !.OR.N.EQ.5

              IF (N.eq.1) A=BMEinlet(NC)+BMR(N,L) 
			if (N.eq.2) A=BMEoutlet(NC)+BMR(N,L)
			if (n.eq.6) a=BMElower(nc)+bmr(N,L)
              !IF (N.GT.2.AND.N.LT.5) A=BME(N,L)*CG(3,NC)/FW+BMR(N,L)
!			IF (N.EQ.5) A=BMELOWER(NC)+BMR(N,L)
              M=A
			BMR(N,L)=A-M
*--M molecules enter, remainder has been reset

              IF (M.GT.0) THEN
                IF (N.EQ.1.OR.N.EQ.2) THEN
                  IF (ABS(vel(1)).GT.1.E-6) THEN !vel(1)
                    IF (N.EQ.1) SC=vel(1)/vmp !VFX/VMP
                    IF (N.EQ.2) SC=-VEL(1)/VMP !-VFX/VMP
                  END IF
                END IF
                IF (n.eq.6) then	!(N.EQ.3.OR.N.EQ.4) THEN
                  IF (ABS(vel(2)).GT.1.E-6) THEN
				  sc=vel(2)/vmp
                    !IF (N.EQ.3) SC=VFY/VMP
                    !IF (N.EQ.4) SC=-VFY/VMP
                  END IF
                END IF
                FS1=SC+SQRT(SC*SC+2.)
                FS2=0.5*(1.+SC*(2.*SC-FS1))
* the above constants are required for the entering distn. of eqn (12.5)
                DO 4 K=1,M
                  IF (NM.LT.MNM) THEN
                    NM=NM+1
*--NM is now the number of the new molecule
                    IF ((N.LT.3.AND.ABS(VEL(1)).GT.1.E-6).OR.
     &                  (N.GT.2.AND.ABS(VEL(2)).GT.1.E-6)) THEN
                      QA=3.
                      IF (SC.LT.-3.) QA=ABS(SC)+1.
2                     U=-QA+2.*QA*RF(0)
*--U is a potential normalised thermal velocity component
                      UN=U+SC
*--UN is a potential inward velocity component
                      IF (UN.LT.0.) GO TO 2
                      A=(2.*UN/FS1)*EXP(FS2-U*U)
                      IF (A.LT.RF(0)) GO TO 2
*--the inward normalised vel. component has been selected (eqn (12.5))
                      IF (N.EQ.1) PV(1,NM)=UN*VMP
                      IF (N.EQ.2) PV(1,NM)=-UN*VMP
                      IF (N.EQ.3) PV(2,NM)=UN*VMP
                      IF (N.EQ.4) PV(2,NM)=-UN*VMP
					if (n.eq.6) pv(2,nm)=un*vmp
                    ELSE
                      IF (N.EQ.1) PV(1,NM)=SQRT(-LOG(RF(0)))*VMP
                      IF (N.EQ.2) PV(1,NM)=-SQRT(-LOG(RF(0)))*VMP
                      IF (N.EQ.3) PV(2,NM)=SQRT(-LOG(RF(0)))*VMP
                      IF (N.EQ.4) PV(2,NM)=-SQRT(-LOG(RF(0)))*VMP
					if (n.eq.6) pv(2,nm)=sqrt(-log(rf(0)))*vmp
*--for a stationary external gas, use eqn (12.3)
                    END IF
                    IF (N.LT.3) THEN
                      CALL RVELC(PV(2,NM),PV(3,NM),VMP)
                      PV(2,NM)=PV(2,NM)+VEL(2)	!VFY
                    END IF
                    IF (N.GT.2) THEN
                      CALL RVELC(PV(1,NM),PV(3,NM),VMP)
                      PV(1,NM)=PV(1,NM)+VEL(1)
                    END IF

*--a single call of RVELC generates the two normal velocity components
                    IF (N.EQ.1) FTMP1=FTMP
				  IF (N.EQ.2) FTMP1=TEMP
				  if (n.eq.6) ftmp1=temp

				  IF (ISPR(1,L).GT.0) CALL SROT(PR(NM),FTMP,ISPR(1,L))
                    IF (N.EQ.1) PP(1,NM)=CB(1)+0.001*CG(3,1)
                    IF (N.EQ.2) PP(1,NM)=cb(8)-0.001*cg(3,1)	!CB(2)-0.001*CG(3,MNC)
                    IF (N.EQ.3) PP(2,NM)=CB(3)+0.0001*(CG(6,1)-CG(5,1))
                    IF (N.EQ.4) PP(2,NM)=CB(4)-0.0001*(CG(6,1)-CG(5,1))
			  if (n.eq.6) pp(2,nm)=cb(9)+0.0001*(cg(6,mnc)-cg(5,mnc))
*--the molecule is moved just off the boundary
                    IPS(NM)=L
                    IF (N.LT.3) THEN
                      IF (N.EQ.1) MC=(NC-1)*NCX+1
                      IF (N.EQ.2) MC=ncx*ncy+nc*nbx	!NC*NCX
					
                      !PP(2,NM)=CG(4,MC)+RF(0)*(CG(6,1)-CG(5,1)) 
					!SINCE IN THE INLET ZONE THE CELLS ARE UNIFORM. If not, as follows
					  XP1=CG(1,MC)
					  YP1=CG(4,MC)
					  XP2=CG(2,MC)
					  YP2=CG(5,MC)
					  MSLOPE=(YP2-YP1)/(XP2-XP1)
					  YL=MSLOPE*(PP(1,NM)-XP1)+YP1

					  XP1=CG(1,MC)
					  YP1=CG(7,MC)
					  XP2=CG(2,MC)
					  YP2=CG(6,MC)
					  MSLOPE=(YP2-YP1)/(XP2-XP1)
					  YU=MSLOPE*(PP(1,NM)-XP1)+YP1

					  CALL Y_LOWER_BOUND(PP(1,NM),YLB)
					  PP(2,NM)=YL+RF(0)*(YU-YL) 
					  CALL Y_LOWER_BOUND(PP(1,NM),YLB)
					  IF(PP(2,NM).LT.YLB)THEN
						PRINT*,'ERROR IN ENTERING MOLECULES'
					  ENDIF
					  if(n.eq.1) then
						NYCH=FLOOR((PP(2,NM)-YLB)/(CB(4)-YLB)*NCY)+1
					  elseif(n.eq.2)then
						NYCH=FLOOR((PP(2,NM)-YLB)/(CB(4)-YLB)*nby)+1
					  endif

					  IF(NYCH.NE.NC)THEN
						PRINT*,'ERROR IN ENTER 2',NC,NYCH
					  ENDIF
                    END IF
                    IF (N.GT.2) THEN
                      IF (N.EQ.3) MC=NC
                      IF (N.EQ.4) MC=(NCY-1)*NCX+NC
					if (n.eq.6) mc=ncx*ncy+nc
                      PP(1,NM)=CG(1,MC)+RF(0)*CG(3,MC)
                    END IF
                    MSCX=((PP(1,NM)-CG(1,MC))/CG(3,MC))*(NSCX-.001)+1
!                    MSCY=((PP(2,NM)-CG(4,MC))/
!     &					(CG(6,MC)-CG(5,MC)))*(NSCY-.001)+1
					!SINCE IN THE INLET ZONE THE CELLS ARE UNIFORM
		CALL ROWLocalSubCell(cg(1,MC),cg(4,MC),cg(2,MC),cg(5,MC),
     &	cg(2,MC),cg(6,MC),cg(1,MC),cg(7,MC),PP(1,NM),PP(1,NM),mscy)
                    MSC=(MSCY-1)*NSCX+MSCX+NSCX*NSCY*(MC-1)
*--MSC is the new sub-cell number
                    IF (MSC.LT.1) MSC=1
                    IF (MSC.GT.MNSC) MSC=MNSC
                    IPL(NM)=MSC
                  ELSE
                    WRITE (*,*) 
     &' WARNING: EXCESS MOLECULE LIMIT - RESTART WITH AN INCREASED FNUM'
                  END IF
4               CONTINUE
              END IF
10          CONTINUE
20        CONTINUE
        END IF
100   CONTINUE
*--now the jet molecules
      IF (IJET.GT.0) THEN
        NCS=LIMJ(3)-LIMJ(2)+1
        DO 150 NC=1,NCS
          DO 120 L=1,MNSP
*--consider each species in turn
            VMP=SQRT(2.*BOLTZ*TMPJ/SP(5,L))
            NCL=NC+LIMJ(2)-1
            IF (IJET.LT.3) A=BMEJ(L)*CG(3,NCL)/WJ+BMRJ(L)
            IF (IJET.GT.2) A=BMEJ(L)*(CG(6,(NCL-1)*NCX+1)-
     &				CG(5,(NCL-1)*NCX+1))/WJ+BMRJ(L)
            M=A
            BMRJ(L)=A-M
*--M molecules enter, remainder has been reset
            IF (M.GT.0) THEN
              IF (ABS(FVJ).GT.1.E-6) SC=FVJ/VMP
              FS1=SC+SQRT(SC*SC+2.)
              FS2=0.5*(1.+SC*(2.*SC-FS1))
* the above constants are required for the entering distn. of eqn (12.5)
              DO 105 K=1,M
                IF (NM.LT.MNM) THEN
                  NM=NM+1
*--NM is now the number of the new molecule
                  IF (ABS(FVJ).GT.1.E-6) THEN
                    QA=3.
                    IF (SC.LT.-3.) QA=ABS(SC)+1.
102                 U=-QA+2.*QA*RF(0)
*--U is a potential normalised thermal velocity component
                    UN=U+SC
*--UN is a potential inward velocity component
                    IF (UN.LT.0.) GO TO 102
                    A=(2.*UN/FS1)*EXP(FS2-U*U)
                    IF (A.LT.RF(0)) GO TO 102
*--the inward normalised vel. component has been selected (eqn (12.5))
                    IF (IJET.EQ.1) PV(2,NM)=UN*VMP
                    IF (IJET.EQ.2) PV(2,NM)=-UN*VMP
                    IF (IJET.EQ.3) PV(1,NM)=UN*VMP
                    IF (IJET.EQ.4) PV(1,NM)=-UN*VMP
                  ELSE
                    IF (IJET.EQ.1) PV(2,NM)=SQRT(-LOG(RF(0)))*VMP
                    IF (IJET.EQ.2) PV(2,NM)=-SQRT(-LOG(RF(0)))*VMP
                    IF (IJET.EQ.3) PV(1,NM)=SQRT(-LOG(RF(0)))*VMP
                    IF (IJET.EQ.4) PV(1,NM)=-SQRT(-LOG(RF(0)))*VMP
*--for a stationary external gas, use eqn (12.3)
                  END IF
                  IF (IJET.LT.3) CALL RVELC(PV(1,NM),PV(3,NM),VMP)
                  IF (IJET.GT.2) CALL RVELC(PV(2,NM),PV(3,NM),VMP)
*--a single call of RVELC generates the two normal velocity components
                  IF (ISPR(1,L).GT.0) CALL SROT(PR(NM),TMPJ,ISPR(1,L))
                  IF (IJET.LT.3) THEN
                    MC=(LIMJ(1)-1)*NCX+LIMJ(2)-1+NC
                    YJ=CG(4,MC)
                    IF (IJET.EQ.2) MC=MC-NCX
                  END IF
                  IF (IJET.GT.2) THEN
                    MC=LIMJ(1)+(LIMJ(2)-1)*NCX+(NC-1)*NCX
                    XJ=CG(1,MC)
                    IF (IJET.EQ.4) MC=MC-1
                  END IF
                  IF (IJET.EQ.1) PP(2,NM)=YJ+0.0001*(CG(6,MC)-CG(5,MC))
                  IF (IJET.EQ.2) PP(2,NM)=YJ-0.0001*(CG(6,MC)-CG(5,MC))
                  IF (IJET.EQ.3) PP(1,NM)=XJ+0.001*CG(3,MC)
                  IF (IJET.EQ.4) PP(1,NM)=XJ-0.001*CG(3,MC)
*--the molecule is moved just off the boundary
                  IPS(NM)=L
                  IF (IJET.LT.3) PP(1,NM)=CG(1,MC)+RF(0)*CG(3,MC)
                  IF (IJET.GT.2) PP(2,NM)=CG(4,MC)+RF(0)
     &								   *(CG(6,MC)-CG(5,MC))
                  MSCX=((PP(1,NM)-CG(1,MC))/CG(3,MC))*(NSCX-.001)+1
                  MSCY=((PP(2,NM)-CG(4,MC))/
     &					(CG(6,MC)-CG(5,MC)))*(NSCY-.001)+1
                  MSC=(MSCY-1)*NSCX+MSCX+NSCX*NSCY*(MC-1)
*--MSC is the new sub-cell number
                  IF (MSC.LT.1) MSC=1
                  IF (MSC.GT.MNSC) MSC=MNSC
                  IPL(NM)=MSC
                ELSE
                  WRITE (*,*) 
     &' WARNING: EXCESS MOLECULE LIMIT - RESTART WITH AN INCREASED FNUM'
                END IF
105           CONTINUE
            END IF
120       CONTINUE
150     CONTINUE
      END IF
      RETURN
      END
!======================================================================
      SUBROUTINE REFLECT2(N,KS,K,XC,YC,MC)
*	 CALL REFLECT2(N,KS,MCS,XC,YSC,MC)
*--reflection of molecule N from surface KS, element K,
*----location XC,YC, cell MC
*
      Include 'common.txt'

	REAL SLOPE,U,V,ANGLE,YLB,M1,ALFA,X1,X2,Y1,Y2,VT,VN
	
	KKK=k
	IF (ISURF(KS)==3) kkk=k+NCX
	!write(*,*) ISURF(KS),kkk
      L=IPS(N)
*--sample the surface properies due to the incident molecules
      CSS(1,K,L)=CSS(1,K,L)+1.
      IF (ISURF(KS).EQ.1) THEN
        CSS(2,K,L)=CSS(2,K,L)-SP(5,L)*PV(2,N)
        CSS(4,K,L)=CSS(4,K,L)+SP(5,L)*PV(1,N)
      END IF
      IF (ISURF(KS).EQ.2) THEN
        CSS(2,K,L)=CSS(2,K,L)+SP(5,L)*PV(2,N)
        CSS(4,K,L)=CSS(4,K,L)+SP(5,L)*PV(1,N)
      END IF
      IF (ISURF(KS).EQ.3) THEN
        CSS(2,K,L)=CSS(2,K,L)-SP(5,L)*PV(1,N)
        CSS(4,K,L)=CSS(4,K,L)+SP(5,L)*PV(2,N)
      END IF
      IF (ISURF(KS).EQ.4) THEN
        CSS(2,K,L)=CSS(2,K,L)+SP(5,L)*PV(1,N)
        CSS(4,K,L)=CSS(4,K,L)+SP(5,L)*PV(2,N)
      END IF
      CSS(5,K,L)=CSS(5,K,L)+0.5*SP(5,L)
     &           *(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
      CSS(7,K,L)=CSS(7,K,L)+PR(N)
*
      IF (TSURF(KS).LT.0.) THEN
*--specular reflection
        IF (ISURF(KS).EQ.1.OR.ISURF(KS).EQ.2)THEN
*	   PV(2,N)=-PV(2,N)	!OLD FORMULA FOR CONSTANT AREA CHANNEL
		IF(XC.LT.XYN(1,1))THEN
			PV(2,N)=-PV(2,N)			
		ELSEIF(XC.GT.XYN(1,1).AND.XC.LT.XYN(2,1))THEN
			X1=XYN(1,1);Y1=XYN(1,2)
			X2=XYN(2,1);Y2=XYN(2,2)
			M1=(Y2-Y1)/(X2-X1)
			ALFA=ATAN(M1)
			VT=PV(1,N)*COS(ALFA)+PV(2,N)*SIN(ALFA)
			VN=-PV(1,N)*SIN(ALFA)+PV(2,N)*COS(ALFA)
			VT=VT
			VN=-VN
			PV(1,N)=VT*COS(ALFA)-VN*SIN(ALFA)
			PV(2,N)=VT*SIN(ALFA)+VN*COS(ALFA)
		ELSEIF(XC.GT.XYN(2,1).AND.XC.LT.XYN(3,1))THEN
			X1=XYN(2,1);Y1=XYN(2,2)
			X2=XYN(3,1);Y2=XYN(3,2)
			M1=(Y2-Y1)/(X2-X1)
			ALFA=ATAN(M1)
			VT=PV(1,N)*COS(ALFA)+PV(2,N)*SIN(ALFA)
			VN=-PV(1,N)*SIN(ALFA)+PV(2,N)*COS(ALFA)
			VT=VT
			VN=-VN
			PV(1,N)=VT*COS(ALFA)-VN*SIN(ALFA)
			PV(2,N)=VT*SIN(ALFA)+VN*COS(ALFA)
	    ELSEIF((XC.GT.XYN(3,1)).and.(xc.lt.cb(2)))THEN
			PV(2,N)=-PV(2,N)
		ENDIF
	  ENDIF
        IF (ISURF(KS).EQ.3.OR.ISURF(KS).EQ.4) PV(1,N)=-PV(1,N)
      ELSE IF (ALPI(KS).LT.0.) THEN
*--diffuse reflection
		if (Tw(KKK).LT.0) then 
			Tw(KKK)=ABS(Tw(KKK))
			write(*,*) 'negative Tw!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
		End if 
        VMP=SQRT(2.*BOLTZ*Tw(KKK)/SP(5,L))
*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
        IF (ISURF(KS).EQ.1) THEN
          V=SQRT(-LOG(RF(0)))*VMP
          CALL RVELC(U,PV(3,N),VMP)

		if(xc.lt.cb(2))then
		IF((XC.LT.XYN(1,1)).OR.(XC.GT.XYN(3,1)))THEN
			PV(1,N)=U
			PV(2,N)=V
		ELSEIF(XC.LT.XYN(2,1).AND.XC.GT.XYN(1,1))THEN
			SLOPE=(XYN(2,2)-XYN(1,2))/(XYN(2,1)-XYN(1,1))
			ANGLE=ATAN(SLOPE)
			PV(1,N)=U*COS(ANGLE)-V*SIN(ANGLE)
			PV(2,N)=U*SIN(ANGLE)+V*COS(ANGLE)
		ELSEIF(XC.LT.XYN(3,1).AND.XC.GT.XYN(2,1))THEN
			SLOPE=(XYN(3,2)-XYN(2,2))/(XYN(3,1)-XYN(2,1))
			ANGLE=ATAN(SLOPE)
			PV(1,N)=U*COS(ANGLE)-V*SIN(ANGLE)
			PV(2,N)=U*SIN(ANGLE)+V*COS(ANGLE)
		ENDIF
		else
			print*,'error in reflect, xc.gt.cb(2)'
		endif
*          PV(2,N)=SQRT(-LOG(RF(0)))*VMP
*         CALL RVELC(PV(1,N),PV(3,N),VMP)
        END IF
        IF (ISURF(KS).EQ.2) THEN
          PV(2,N)=-SQRT(-LOG(RF(0)))*VMP
          CALL RVELC(PV(1,N),PV(3,N),VMP)
        END IF
        IF (ISURF(KS).EQ.3) THEN
          PV(1,N)=SQRT(-LOG(RF(0)))*VMP
          CALL RVELC(PV(2,N),PV(3,N),VMP)
        END IF
        IF (ISURF(KS).EQ.4) THEN
          PV(1,N)=-SQRT(-LOG(RF(0)))*VMP
          CALL RVELC(PV(2,N),PV(3,N),VMP)
        END IF
*--the normal velocity component has been generated (eqn(12.3))
*--a single call of RVELC generates the two tangential vel. components
        IF (ISPR(1,L).GT.0) CALL SROT(PR(N),Tw(KKK),ISPR(1,L))
      ELSE IF (ALPI(KS).GE.0) THEN
*--Cercignani-Lampis-Lord reflection model
        VMP=SQRT(2.*BOLTZ*Tw(KKK)/SP(5,L))
*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
        IF (ISURF(KS).EQ.1.OR.ISURF(KS).EQ.2) THEN
          IF (ISURF(KS).EQ.1) VNI=-PV(2,N)/VMP
          IF (ISURF(KS).EQ.2) VNI=PV(2,N)/VMP
          UPI=PV(1,N)/VMP
        END IF
        IF (ISURF(KS).EQ.3.OR.ISURF(KS).EQ.4) THEN
          IF (ISURF(KS).EQ.3) VNI=-PV(1,N)/VMP
          IF (ISURF(KS).EQ.4) VNI=PV(1,N)/VMP
          UPI=PV(2,N)/VMP
        END IF
        WPI=PV(3,N)/VMP
        ANG=ATAN2(WPI,UPI)
        VPI=SQRT(UPI*UPI+WPI*WPI)
*--VNI is the normalized incident normal vel. component (always +ve)
*--VPI is the normalized incident tangential vel. comp. in int. plane
*--ANG is the angle between the interaction plane and the x or y axis
*--first the normal component
        ALPHAN=ALPN(KS)
        R=SQRT(-ALPHAN*LOG(RF(0)))
        TH=2.*PI*RF(0)
        UM=SQRT(1.-ALPHAN)*VNI
        VN=SQRT(R*R+UM*UM+2.*R*UM*COS(TH))
*--VN is the normalized magnitude of the reflected normal vel. comp. from eqns (14.3)
*--then the tangential component
        ALPHAT=ALPT(KS)*(2.-ALPT(KS))
        R=SQRT(-ALPHAT*LOG(RF(0)))
        TH=2.*PI*RF(0)
        UM=SQRT(1.-ALPHAT)*VPI
        VP=UM+R*COS(TH)
        WP=R*SIN(TH)
*--VP,WP are the normalized reflected tangential vel. components in and
*----normal to the interaction plane, from eqns (14.4) and (14.5)
        IF (ISURF(KS).EQ.1.OR.ISURF(KS).EQ.2) THEN
          IF (ISURF(KS).EQ.1) PV(2,N)=VN*VMP
          IF (ISURF(KS).EQ.2) PV(2,N)=-VN*VMP
          PV(1,N)=(VP*COS(ANG)-WP*SIN(ANG))*VMP
        END IF
        IF (ISURF(KS).EQ.3.OR.ISURF(KS).EQ.4) THEN
          IF (ISURF(KS).EQ.3) PV(1,N)=VN*VMP
          IF (ISURF(KS).EQ.4) PV(1,N)=-VN*VMP
          PV(2,N)=(VP*COS(ANG)-WP*SIN(ANG))*VMP
        END IF
        PV(3,N)=(VP*SIN(ANG)+WP*COS(ANG))*VMP
        IF (ISPR(1,L).GT.0) THEN
*--set CLL rotational energy by analogy with normal vel. component
          ALPHAI=ALPI(KS)
          OM=SQRT(PR(N)*(1.-ALPHAI)/(BOLTZ*Tw(KKK)))
          IF (ISPR(1,L).EQ.2) THEN
            R=SQRT(-ALPHAI*LOG(RF(0)))
            CTH=COS(2.*PI*RF(0))
          ELSE
*--for polyatomic case, apply acceptance-rejection based on eqn (14.6)
10          X=4.*RF(0)
            A=2.7182818*X*X*EXP(-X*X)
            IF (A.LT.RF(0)) GO TO 10
            R=SQRT(ALPHAI)*X
            CTH=2.*RF(0)-1.
          END IF
          PR(N)=BOLTZ*Tw(KKK)*(R*R+OM*OM+2.*R*OM*CTH)
        END IF
      END IF
      
	IF (ISURF(KS).EQ.1) THEN
        PP(1,N)=XC
        PP(2,N)=YC+0.001*(CG(6,MC)-CG(5,MC))
	  CALL Y_LOWER_BOUND(XC,YLB)
	  IF(PP(2,N).LT.YLB) PRINT*,'**ERROR IN REFLECT==>PP(2) < YLB)'
      END IF
      IF (ISURF(KS).EQ.2) THEN
        PP(1,N)=XC
        PP(2,N)=YC-0.001*(CG(6,MC)-CG(5,MC))
      END IF
      IF (ISURF(KS).EQ.3) THEN
        PP(1,N)=XC+0.001*CG(3,MC)
        PP(2,N)=YC
      END IF
      IF (ISURF(KS).EQ.4) THEN
        PP(1,N)=XC-0.001*CG(3,MC)
        PP(2,N)=YC
      END IF
      IPL(N)=(MC-1)*NSCX*NSCY+1
*--sample the surface properties due to the reflected molecules
      IF (ISURF(KS).EQ.1) CSS(3,K,L)=CSS(3,K,L)+SP(5,L)*PV(2,N)
      IF (ISURF(KS).EQ.2) CSS(3,K,L)=CSS(3,K,L)-SP(5,L)*PV(2,N)
      IF (ISURF(KS).EQ.3) CSS(3,K,L)=CSS(3,K,L)+SP(5,L)*PV(1,N)
      IF (ISURF(KS).EQ.4) CSS(3,K,L)=CSS(3,K,L)-SP(5,L)*PV(1,N)
      IF (ISURF(KS).EQ.1) CSS(9,K,L)=CSS(9,K,L)-SP(5,L)*PV(1,N)
      IF (ISURF(KS).EQ.2) CSS(9,K,L)=CSS(9,K,L)-SP(5,L)*PV(1,N)
      IF (ISURF(KS).EQ.3) CSS(9,K,L)=CSS(9,K,L)-SP(5,L)*PV(2,N)
      IF (ISURF(KS).EQ.4) CSS(9,K,L)=CSS(9,K,L)-SP(5,L)*PV(2,N)
      CSS(6,K,L)=CSS(6,K,L)-0.5*SP(5,L)
     &           *(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
      CSS(8,K,L)=CSS(8,K,L)-PR(N)
      RETURN
      END
!======================================================================
      SUBROUTINE REMOVE(N)
*
*--remove molecule N and replace it by molecule NM
*
	Include 'common.txt'

      PP(1,N)=PP(1,NM)
      PP(2,N)=PP(2,NM)
      DO 100 M=1,3
        PV(M,N)=PV(M,NM)
100   CONTINUE
      PR(N)=PR(NM)
      IPL(N)=IPL(NM)
      IPS(N)=IPS(NM)
      NM=NM-1
      N=N-1
      RETURN
      END
!======================================================================
      SUBROUTINE SAMPI2
*
*--initialises all the sampling variables
*
	Include 'common.txt'

      NSMP=0
      TIMI=TIME
      DO 200 L=1,MNSP
        DO 50 N=1,MNC
          CS(1,N,L)=1.E-6
          DO 20 M=2,7
            CS(M,N,L)=0.
20        CONTINUE
          CSR(N,L)=0.
50      CONTINUE
        DO 100 N=1,MNSE
          CSS(1,N,L)=1.E-6
          DO 60 M=2,9
            CSS(M,N,L)=0.
60        CONTINUE
100     CONTINUE
200   CONTINUE
      RETURN
      END
!======================================================================
      SUBROUTINE SAMPLE2
*
*--sample the molecules in the flow.
*
	Include 'common.txt'

      NSMP=NSMP+1
      DO 100 NN=1,MNSG
        DO 50 N=1,MNC
          L=IC(2,N,NN)
          IF (L.GT.0) THEN
            DO 10 J=1,L
              K=IC(1,N,NN)+J
              M=IR(K)
              I=IPS(M)
              CS(1,N,I)=CS(1,N,I)+1
              DO 5 LL=1,3
                CS(LL+1,N,I)=CS(LL+1,N,I)+PV(LL,M)
                CS(LL+4,N,I)=CS(LL+4,N,I)+PV(LL,M)**2
5             CONTINUE
              CSR(N,I)=CSR(N,I)+PR(M)
10          CONTINUE
          END IF
50      CONTINUE
100   CONTINUE
      RETURN
      END
!======================================================================
      SUBROUTINE OUT2
*
*--output a progressive set of results to file DSMC2.OUT.
*
      Include 'common.txt'
	INCLUDE 'PROPERTY.TXT'

	real,allocatable::attribute(:,:)
	real::prop(10),dx,dye,Ch1,Ch2,Ch3,Cy1,Cy2,Rm1,Rm2,Rm3,Rm4,Rm5
     &,Sumkn,sumt
	integer ::i, Mc,Ncell

	allocate(attribute(mnc,10))

	DBOLTZ=BOLTZ
	
	QWALL=0
	TAWALL=0	

	KS=1
	L=1
     	
	SumKn=0.
	SumT=0.

	DO 400 N=1,MNC
		CALL PROPERTIES (N)
		attribute(n,1)=den
		attribute(n,2)=tt
		attribute(n,3)=trot
		attribute(n,4)=temp
		attribute(n,5)=vel(1)
		If (n.LT.NCX+1) then 
			 Vx(n)=vel(1)
		End if 
		attribute(n,6)=vel(2)
		attribute(n,7)=vel(3)
		attribute(n,8)=xm
		attribute(n,9)=p
		attribute(n,10)=Kn
400	END DO

	!Mass Flow rate and Thrust 
	Ncell=NCX*NCY
	do n=1,Ncell
		SumKn=SumKn+attribute(n,10)
		Sumt=sumt+attribute(n,4)
	end do 
	Sumt=sumt/Ncell
	Sumkn=Sumkn/ncell
	Tavg=Sumt

	Rm1=0.
	Rm2=0.
	Rm3=0.
	Rm4=0.
	Rm5=0.

	
	
	DO N=2,4
		
		DO NC=1,NCY
			IF (n.EQ.2) Mc=(NC-1)*Ncx+0.15*Ncx
			IF (n.EQ.3) Mc=(NC-1)*Ncx+0.25*Ncx
			IF (n.EQ.4) Mc=(NC-1)*Ncx+0.60*Ncx
			
			IF (N==2) THEN
				Rm2=Rm2+attribute(Mc,1)*attribute(Mc,5)
			ELSEIF (N==3) THEN
				Rm3=Rm3+attribute(Mc,1)*attribute(Mc,5)
			ELSEIF (N==4) THEN
				Rm4=Rm4+attribute(Mc,1)*attribute(Mc,5)
			ENDIF
		ENDDO
	ENDDO
	
	dx=(cb(2)-cb(1))/NCX	
	dye=(cb(4)-cb(3))/NCY

	

	Ch3=CB(4)	!outlet height
	Ch2=CB(4)-CB(5) !Throat height
	Ch1=CB(6)	!inlet height
	Cy1=(Ch1-Ch2)/(INZ(2)*dx)*(0.25-0.15)*NCX*dx+Ch2
	Cy2=(Ch3-Ch2)/((INZ(2)-INZ(3))*dx)*(0.25-0.60)*NCX*dx+Ch2

	!open(123,file='test.txt')
	!write(123,*) Ch2
	!close (123)

	Rm2=Rm2/NCY*2.*22.6e-6
	Rm3=Rm3/NCY*2.*Ch2
	Rm4=Rm4/NCY*2.*40.2e-6
	Rmdot=Rm3
	Rm5=Rm4
	!thrust calculation 
	Do i=1,NCY
		Mc=i*NCX
	Rm1=Rm1+2*Abs(attribute(Mc,1)*attribute(Mc,5))*attribute(Mc,5)*dye
     &+2*attribute(Mc,9)*dye
	end do
	
	
	OPEN (30,FILE='Convergence.plt',position="append")
		WRITE (30,99000) NPR,Rm1,Rm2,Rm3,Rm4,Rm5,sumt,sumkn	
	CLOSE (30)

99000	FORMAT (I8,7ES20.4)
*	
	!End mass flow rate and thrust 
	
	OPEN (9,FILE='Contour_2_zone.PLT')
	WRITE(9,*)'TITLE = "Subsonic"'
      WRITE (9,*)'Variables=X,Y,Density, TrTemp, RotTemp, T,U,V,W
     &,Mach,Pressure,Kundsen'
	WRITE(9,*)'ZONE T="MAIN ZONE", I=',ncx+1,', J=',NCY+1,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	DO j=1,ncy+1
	  do i=1,ncx+1
		if(i.eq.1)then
			if(j.eq.1)then
				prop(:)=attribute(1,:)
			elseif(j.eq.(ncy+1))then
				prop(:)=attribute((ncy-1)*ncx+1,:)
			else
				prop(:)=(attribute((j-2)*ncx+1,:) +
     &				     attribute((j-1)*ncx+1,:))/2.
			endif
		elseif(i.eq.(ncx+1))then
			if(j.eq.1)then
				prop(:)=attribute(ncy,:)
			elseif(j.eq.(ncy+1))then
				prop(:)=attribute(ncy*ncx,:)
			else
				prop(:)=(attribute((j-1)*ncx,:) +
     &				     attribute(   j *ncx,:))/2.
			endif
		else
			if(j.eq.1)then
				prop(:)=(attribute(i-1,:) + 
     &					 attribute( i ,:))/2.
			elseif(j.eq.(ncy+1))then
				prop(:)=(attribute((ncy-1)*ncx+i-1,:) +
     &					 attribute((ncy-1)*ncx+ i ,:))/2.
			else
				prop(:)=(attribute((j-2)*ncx+i-1,:) +
     &					 attribute((j-2)*ncx+i  ,:) +
     &				     attribute((j-1)*ncx+i-1,:) +
     &				     attribute((j-1)*ncx+i  ,:))/4.
			endif
		endif

	    WRITE (9,'(12ES20.4)') xgrid(i,j),ygrid(i,j),prop
!		DEN,TT,TROT,TEMP,VEL(1),VEL(2),VEL(3)
!     &,XM,P/PIN,Kn

	  enddo
	END DO
      
	WRITE(9,*)'ZONE T="BUFFER ZONE", I=',nbx+1,', J=',NbY+1,'
     &,ZONETYPE=Ordered, DATAPACKING=POINT'

	mncm=ncx*ncy
	DO j=1,nby+1
	  do i=1,nbx+1
		if(i.eq.1)then
			if(j.eq.1)then
				prop(:)=attribute(mncm+1,:)
			elseif(j.eq.(nby+1))then
				prop(:)=attribute(mncm+(nby-1)*nbx+1,:)
			else
				prop(:)=(attribute(mncm+(j-2)*nbx+1,:) +
     &				     attribute(mncm+(j-1)*nbx+1,:))/2.
			endif
		elseif(i.eq.(nbx+1))then
			if(j.eq.1)then
				prop(:)=attribute(mncm+nby,:)
			elseif(j.eq.(nby+1))then
				prop(:)=attribute(mncm+nby*nbx,:)
			else
				prop(:)=(attribute(mncm+(j-1)*nbx,:) +
     &				     attribute(mncm+   j *nbx,:))/2.
			endif
		else
			if(j.eq.1)then
				prop(:)=(attribute(mncm+i-1,:) + 
     &					 attribute(mncm+ i ,:))/2.
			elseif(j.eq.(nby+1))then
				prop(:)=(attribute(mncm+(nby-1)*nbx+i-1,:) +
     &					 attribute(mncm+(nby-1)*nbx+ i ,:))/2.
			else
				prop(:)=(attribute(mncm+(j-2)*nbx+i-1,:) +
     &					 attribute(mncm+(j-2)*nbx+i  ,:) +
     &					 attribute(mncm+(j-1)*nbx+i-1,:) +
     &					 attribute(mncm+(j-1)*nbx+i  ,:))/4.
			endif
		endif

	    WRITE (9,'(12ES20.4)') x_b(i,j),y_b(i,j),prop
!		DEN,TT,TROT,TEMP,VEL(1),VEL(2),VEL(3)
!     &,XM,P/PIN,Kn

	  enddo
	END DO
	
	CLOSE (9)

	END SUBROUTINE
!======================================================================
      SUBROUTINE PROPERTIES (I)

      Include 'common.txt'
	INCLUDE 'PROPERTY.TXT'


	N=I
	
	DBOLTZ=BOLTZ

!	DO 400 N=1,MNC
        A=FNUM/(CC(N)*NSMP)
        SN=0.
        SM=0.
        DO 250 K=1,3
          SMU(K)=0.
250     CONTINUE
        SMCC=0.
        SRE=0.
        SRDF=0.
        DO 300 L=1,MNSP
          SN=SN+CS(1,N,L)
*--SN is the number sum
          SM=SM+SP(5,L)*CS(1,N,L)
*--SM is the sum of molecular masses
          DO 260 K=1,3
            SMU(K)=SMU(K)+SP(5,L)*CS(K+1,N,L)
*--SMU(1 to 3) are the sum of mu, mv, mw
260       CONTINUE
          SMCC=SMCC+(CS(5,N,L)+CS(6,N,L)+CS(7,N,L))*SP(5,L)
*--SMCC is the sum of m(u**2+v**2+w**2)
          SRE=SRE+CSR(N,L)
*--SRE is the sum of rotational energy
          SRDF=SRDF+ISPR(1,L)*CS(1,N,L)
*--SRDF is the sum of the rotational degrees of freedom
          SUU=SUU+SP(5,L)*CS(5,N,L)
*--SUU is the sum of m*u*u
300     CONTINUE
        DENN=SN*A
*--DENN is the number density, see eqn (1.34)
	 
	  mfp=1.0/(SQRT(2.0)*(pi*sp(1,1)**2)*DENN)  
        kn=mfp/CB(4)

        DEN=DENN*SM/SN
*--DEN is the density, see eqn (1.42)
        DO 350 K=1,3
          VEL(K)=SMU(K)/SM
          SVEL(K,N)=VEL(K)
350     CONTINUE
*--VEL and SVEL are the stream velocity components, see eqn (1.43)
        UU=VEL(1)**2+VEL(2)**2+VEL(3)**2
        TT=(SMCC-SM*UU)/(3.D00*DBOLTZ*SN)
*--TT is the translational temperature, see eqn (1.51)
        IF (SRDF.GT.1.E-6) TROT=(2.D00/DBOLTZ)*SRE/SRDF
*--TROT is the rotational temperature, see eqn (11.11)
        TEMP=(3.D00*TT+(SRDF/SN)*TROT)/(3.+SRDF/SN)
	  IF (TEMP.EQ.0) TEMP=300.
	  P=DENN*BOLTZ*TEMP
*--TEMP is the overall temperature, see eqn (11.12)
        CT(N)=TEMP
        XC=0.5*(CG(1,N)+CG(2,N))
        YC=0.5*(CG(4,N)+CG(5,N))
	  GAMA=(5.+ISPR(1,1))/(3.+ISPR(1,1))
	  GASR=BOLTZ/SP(5,1)
	  ASOUND=SQRT(GAMA*GASR*TEMP)
	  XM=SQRT(UU)/ASOUND
*--XC,YC are the x,y coordinates of the midpoint of the cell

	! Corrections FOR OUTPUT CELLS
	 if(n.gt.(ncx*ncy))then
	 mc=n-ncx*ncy
	 IF ((mc-(mc/nbx)*nbx).EQ.0) THEN

		 COR_U=(P-POUT)/(DEN*ASOUND)
		 VEL(1)=VEL(1)+COR_U
		 DEN_COR=DEN+(POUT-P)/ASOUND**2
		 T_COR=POUT/(DEN_COR*GASR)
		 DEN=DEN_COR
	     TEMP=T_COR
	 
		 DENN=DEN/SP(5,1)

!		 P=POUT
		 !P=DENN*BOLTZ*TEMP
	 
	 END IF
	 if(mc.lt.nbx)then

		 COR_v=(P-POUT)/(DEN*ASOUND)
		 VEL(2)=VEL(2)+COR_v
		 DEN_COR=DEN+(POUT-P)/ASOUND**2
		 T_COR=POUT/(DEN_COR*GASR)
		 DEN=DEN_COR
	     TEMP=T_COR
	 
		 DENN=DEN/SP(5,1)

!		 P=POUT
		
	 endif

	! Corrections FOR INPUT CELLS
	 elseif(n.lt.(ncx*ncy))then
	 IF ((N-(N/NCX)*NCX).EQ.1) THEN

		COR_U2=(PIN-P)/(DEN*ASOUND)
		VEL(1)=VEL(1)+COR_U2
		
		TEMP=FTMP
!		P=PIN
		DEN=PIN/(GASR*FTMP)
		DENN=DEN/SP(5,1)
		
	!	P=DENN*BOLTZ*TEMP		 		 
	 
	 END IF
	endif

      RETURN
      END SUBROUTINE
!======================================================================
      SUBROUTINE SROT(PR,TEMP,IDF)
*--selects a typical equuilibrium value of the rotational energy PR at
*----the temperature TEMP in a gas with IDF rotl. deg. of f.
*
      COMMON /CONST / PI,SPI,BOLTZ
 
      IF (IDF.EQ.2) THEN
        PR=-LOG(RF(0))*BOLTZ*TEMP
*--for 2 degrees of freedom, the sampling is directly from eqn (11.22)
      ELSE
*--otherwise apply the acceptance-rejection method to eqn (11.23)
        A=0.5*IDF-1.
50      ERM=RF(0)*10.
*--the cut-off internal energy is 10 kT
        B=((ERM/A)**A)*EXP(A-ERM)
        IF (B.LT.RF(0)) GO TO 50
        PR=ERM*BOLTZ*TEMP
      END IF
      RETURN
      END
!======================================================================
      FUNCTION ERF(S)
*
*--calculates the error function of S
*
      B=ABS(S)
      IF (B.GT.4.) THEN
        D=1.
      ELSE
        C=EXP(-B*B)
        T=1./(1.+0.3275911*B)
        D=1.-(0.254829592*T-0.284496736*T*T+1.421413741*T*T*T-
     &    1.453152027*T*T*T*T+1.061405429*T*T*T*T*T)*C
      END IF
      IF (S.LT.0.) D=-D
      ERF=D
      RETURN
      END
!======================================================================
      SUBROUTINE INDEXM
*
*--the NM molecule numbers are arranged in order of the molecule groups
*--and, within the groups, in order of the cells and, within the cells,
*--in order of the sub-cells
*
      Include 'common.txt'
*
      DO 200 MM=1,MNSG
        IG(2,MM)=0
        DO 50 NN=1,MNC
          IC(2,NN,MM)=0
50      CONTINUE
        DO 100 NN=1,MNSC
          ISCG(2,NN,MM)=0
100     CONTINUE
200   CONTINUE
      DO 300 N=1,NM
        LS=IPS(N)
        MG=ISP(LS)
        IG(2,MG)=IG(2,MG)+1
        MSC=IPL(N)
        ISCG(2,MSC,MG)=ISCG(2,MSC,MG)+1
        MC=ISC(MSC)
        IC(2,MC,MG)=IC(2,MC,MG)+1
300   CONTINUE
*--number in molecule groups in the cells and sub-cells have been counte
      M=0
      DO 400 L=1,MNSG
        IG(1,L)=M
*--the (start address -1) has been set for the groups
        M=M+IG(2,L)
400   CONTINUE
      DO 600 L=1,MNSG
        M=IG(1,L)
        DO 450 N=1,MNC
          IC(1,N,L)=M
          M=M+IC(2,N,L)
450     CONTINUE
*--the (start address -1) has been set for the cells
        M=IG(1,L)
        DO 500 N=1,MNSC
          ISCG(1,N,L)=M
          M=M+ISCG(2,N,L)
          ISCG(2,N,L)=0
500     CONTINUE
600   CONTINUE
*--the (start address -1) has been set for the sub-cells
 
      DO 700 N=1,NM
        LS=IPS(N)
        MG=ISP(LS)
        MSC=IPL(N)
        ISCG(2,MSC,MG)=ISCG(2,MSC,MG)+1
        K=ISCG(1,MSC,MG)+ISCG(2,MSC,MG)
        IR(K)=N
*--the molecule number N has been set in the cross-reference array
700   CONTINUE
      RETURN
      END
!======================================================================
      SUBROUTINE SELECT
*--selects a potential collision pair and calculates the product of the
*--collision cross-section and relative speed
*
      Include 'common.txt'
      COMMON /ELAST / VRC(3),VRR,VR,L,M,LS,MS,CVR,MM,NN,N
*
      K=INT(RF(0)*(IC(2,N,NN)-0.001))+IC(1,N,NN)+1
      L=IR(K)
*--the first molecule L has been chosen at random from group NN in cell
100   MSC=IPL(L)
      IF ((NN.EQ.MM.AND.ISCG(2,MSC,MM).EQ.1).OR.
     &    (NN.NE.MM.AND.ISCG(2,MSC,MM).EQ.0)) THEN
*--if MSC has no type MM molecule find the nearest sub-cell with one
        NST=1
        NSG=1
150     INC=NSG*NST
        NSG=-NSG
        NST=NST+1
        MSC=MSC+INC
        IF (MSC.LT.1.OR.MSC.GT.MNSC) GO TO 150
        IF (ISC(MSC).NE.N.OR.ISCG(2,MSC,MM).LT.1) GO TO 150
      END IF
*--the second molecule M is now chosen at random from the group MM
*--molecules that are in the sub-cell MSC
      K=INT(RF(0)*(ISCG(2,MSC,MM)-0.001))+ISCG(1,MSC,MM)+1
      M=IR(K)
      IF (L.EQ.M) GO TO 100
*--choose a new second molecule if the first is again chosen
*
      DO 200 K=1,3
        VRC(K)=PV(K,L)-PV(K,M)
200   CONTINUE
*--VRC(1 to 3) are the components of the relative velocity
      VRR=VRC(1)**2+VRC(2)**2+VRC(3)**2
      VR=SQRT(VRR)
*--VR is the relative speed
      LS=IPS(L)
      MS=IPS(M)
      CVR=VR*SPM(1,LS,MS)*((2.*BOLTZ*SPM(2,LS,MS)/(SPM(5,LS,MS)*VRR))
     &    **(SPM(3,LS,MS)-0.5))/SPM(6,LS,MS)
*--the collision cross-section is based on eqn (4.63)
      RETURN
      END
!======================================================================
      SUBROUTINE ELASTIC
*
*--generate the post-collision velocity components.
*
      Include 'common.txt'
      COMMON /ELAST / VRC(3),VRR,VR,L,M,LS,MS,CVR,MM,NN,N
*
      DIMENSION VRCP(3),VCCM(3)
*--VRCP(3) are the post-collision components of the relative velocity
*--VCCM(3) are the components of the centre of mass velocity
*
      RML=SPM(5,LS,MS)/SP(5,MS)
      RMM=SPM(5,LS,MS)/SP(5,LS)
      DO 100 K=1,3
        VCCM(K)=RML*PV(K,L)+RMM*PV(K,M)
100   CONTINUE
*--VCCM defines the components of the centre-of-mass velocity, eqn (2.1)
      IF (ABS(SPM(4,LS,MS)-1.).LT.1.E-3) THEN
*--use the VHS logic
        B=2.*RF(0)-1.
*--B is the cosine of a random elevation angle
        A=SQRT(1.-B*B)
        VRCP(1)=B*VR
        C=2.*PI*RF(0)
*--C is a random azimuth angle
        VRCP(2)=A*COS(C)*VR
        VRCP(3)=A*SIN(C)*VR
      ELSE
*--use the VSS logic
        B=2.*(RF(0)**SPM(4,LS,MS))-1.
*--B is the cosine of the deflection angle for the VSS model, eqn (11.8)
        A=SQRT(1.-B*B)
        C=2.*PI*RF(0)
        OC=COS(C)
        SC=SIN(C)
        D=SQRT(VRC(2)**2+VRC(3)**2)
        IF (D.GT.1.E-6) THEN
          VRCP(1)=B*VRC(1)+A*SC*D
          VRCP(2)=B*VRC(2)+A*(VR*VRC(3)*OC-VRC(1)*VRC(2)*SC)/D
          VRCP(3)=B*VRC(3)-A*(VR*VRC(2)*OC+VRC(1)*VRC(3)*SC)/D
        ELSE
          VRCP(1)=B*VRC(1)
          VRCP(2)=A*OC*VRC(1)
          VRCP(3)=A*SC*VRC(1)
        END IF
*--the post-collision rel. velocity components are based on eqn (2.22)
      END IF
*--VRCP(1 to 3) are the components of the post-collision relative vel.
      DO 200 K=1,3
        PV(K,L)=VCCM(K)+VRCP(K)*RMM
        PV(K,M)=VCCM(K)-VRCP(K)*RML
200   CONTINUE
      RETURN
      END
!======================================================================
      SUBROUTINE RVELC(U,V,VMP)

*--generates two random velocity components U an V in an equilibrium
*--gas with most probable speed VMP  (based on eqns (C10) and (C12))

      A=SQRT(-LOG(RF(0)))
      B=6.283185308*RF(0)
      U=A*SIN(B)*VMP
      V=A*COS(B)*VMP
      RETURN
      END
!======================================================================
*   GAM.FOR

      FUNCTION GAM(X)

*--calculates the Gamma function of X.

      A=1.
      Y=X
      IF (Y.LT.1.) THEN
        A=A/Y
      ELSE
50      Y=Y-1
        IF (Y.GE.1.) THEN
          A=A*Y
          GO TO 50
        END IF
      END IF
      GAM=A*(1.-0.5748646*Y+0.9512363*Y**2-0.6998588*Y**3+
     &    0.4245549*Y**4-0.1010678*Y**5)
      RETURN
      END
!======================================================================
      SUBROUTINE COLLMR

*--calculates collisions appropriate to DTM in a gas mixture

      Include 'common.txt'
      COMMON /ELAST / VRC(3),VRR,VR,L,M,LS,MS,CVR,MM,NN,N
*
*--VRC(3) are the pre-collision components of the relative velocity
*
      DO 100 N=1,MNC
*--consider collisions in cell N
        DO 50 NN=1,MNSG
          DO 20 MM=1,MNSG
            SN=0.
            DO 10 K=1,MNSP
              IF (ISP(K).EQ.MM) SN=SN+CS(1,N,K)
10          CONTINUE
            IF (SN.GT.1.) THEN
              AVN=SN/FLOAT(NSMP)
            ELSE
              AVN=IC(2,N,MM)
            END IF
*--AVN is the average number of group MM molecules in the cell
            ASEL=0.5*IC(2,N,NN)*AVN*FNUM*CCG(1,N,NN,MM)*DTM/CC(N)
     &           +CCG(2,N,NN,MM)
*--ASEL is the number of pairs to be selected, see eqn (11.5)
            NSEL=ASEL
            CCG(2,N,NN,MM)=ASEL-NSEL
            IF (NSEL.GT.0) THEN
              IF (((NN.NE.MM).AND.(IC(2,N,NN).LT.1.OR.IC(2,N,MM).LT.1))
     &            .OR.((NN.EQ.MM).AND.(IC(2,N,NN).LT.2))) THEN
                CCG(2,N,NN,MM)=CCG(2,N,NN,MM)+NSEL
*--if there are insufficient molecules to calculate collisions,
*--the number NSEL is added to the remainer CCG(2,N,NN,MM)
              ELSE
                CVM=CCG(1,N,NN,MM)
                SELT=SELT+NSEL
                DO 12 ISEL=1,NSEL
*
                  CALL SELECT
*
                  IF (CVR.GT.CVM) CVM=CVR
*--if necessary, the maximum product in CVM is upgraded
                  IF (RF(0).LT.CVR/CCG(1,N,NN,MM)) THEN
*--the collision is accepted with the probability of eqn (11.6)
                    NCOL=NCOL+1
                    SEPT=SEPT+
     &                   SQRT((PP(1,L)-PP(1,M))**2+(PP(2,L)-PP(2,M))**2)
                    COL(LS,MS)=COL(LS,MS)+1.D00
                    COL(MS,LS)=COL(MS,LS)+1.D00
*
                    IF (ISPR(1,LS).GT.0.OR.ISPR(1,MS).GT.0) CALL INELR
*--bypass rotational redistribution if both molecules are monatomic
*
                    CALL ELASTIC
*
                  END IF
12              CONTINUE
                CCG(1,N,NN,MM)=CVM
              END IF
            END IF
20        CONTINUE
50      CONTINUE
100   CONTINUE
 
      RETURN
      END
!======================================================================
      SUBROUTINE INELR
*
*--adjustment of rotational energy in a collision
*
      PARAMETER (MNM=1000000,MNC=3000,MNSC=12000,MNSP=1,MNSG=1,MNSE=90)
*
      COMMON /MOLSR / PR(MNM)
      COMMON /GAS   / SP(5,MNSP),SPM(6,MNSP,MNSP),ISP(MNSP)
      COMMON /GASR  / SPR(3,MNSP,MNSP),ISPR(3,MNSP),CT(MNC)
      COMMON /ELAST / VRC(3),VRR,VR,L,M,LS,MS,CVR,MM,NN,N
*
      DIMENSION IR(2)
*--IR is the indicator for the rotational redistribution
      ETI=0.5*SPM(5,LS,MS)*VRR
*--ETI is the initial translational energy
      ECI=0.
*--ECI is the initial energy in the active rotational modes
      ECF=0.
*--ECF is the final energy in these modes
      ECC=ETI
*--ECC is the energy to be divided
      XIB=2.5-SPM(3,LS,MS)
*--XIB is th number of modes in the redistribution
      IRT=0
*--IRT is 0,1 if no,any redistribution is made
      DO 100 NSP=1,2
*--consider the molecules in turn
        IF (NSP.EQ.1) THEN
          K=L
          KS=LS
          JS=MS
        ELSE
          K=M
          KS=MS
          JS=LS
        END IF
        IR(NSP)=0
        IF (ISPR(1,KS).GT.0) THEN
          IF (ISPR(2,KS).EQ.0) THEN
            ATK=1./SPR(1,KS,JS)
          ELSE
            ATK=1./(SPR(1,KS,JS)+SPR(2,KS,JS)*CT(N)+SPR(3,KS,JS)*CT(N)
     &          **2)
          END IF
*--ATK is the probability that rotation is redistributed to molecule L
          IF (ATK.GT.RF(0)) THEN
            IRT=1
            IR(NSP)=1
            ECC=ECC+PR(K)
            ECI=ECI+PR(K)
            XIB=XIB+0.5*ISPR(1,KS)
          END IF
        END IF
100   CONTINUE
*--apply the general Larsen-Borgnakke distribution function
      IF (IRT.EQ.1) THEN
        DO 150 NSP=1,2
          IF (IR(NSP).EQ.1) THEN
            IF (NSP.EQ.1) THEN
              K=L
              KS=LS
            ELSE
              K=M
              KS=MS
            END IF
            XIB=XIB-0.5*ISPR(1,KS)
*--the current molecule is removed from the total modes
            IF (ISPR(1,KS).EQ.2) THEN
              ERM=1.-RF(0)**(1./XIB)
            ELSE
              XIA=0.5*ISPR(1,KS)
              CALL LBS(XIA-1.,XIB-1.,ERM)
            END IF
            PR(K)=ERM*ECC
            ECC=ECC-PR(K)
*--the available energy is reduced accordingly
            ECF=ECF+PR(K)
          END IF
150     CONTINUE
        ETF=ETI+ECI-ECF
*--ETF  is the post-collision translational energy
*--adjust VR and, for the VSS model, VRC for the change in energy
        A=SQRT(2.*ETF/SPM(5,LS,MS))
        IF (ABS(SPM(4,LS,MS)-1.).LT.1.E-3) THEN
          VR=A
        ELSE
          DO 160 K=1,3
            VRC(K)=VRC(K)*A/VR
160       CONTINUE
          VR=A
        END IF
      END IF
      RETURN
      END
!======================================================================
      SUBROUTINE LBS(XMA,XMB,ERM)
*--selects a Larsen-Borgnakke energy ratio using eqn (11.9)
100   ERM=RF(0)
      IF (XMA.LT.1.E-6.OR.XMB.LT.1.E-6) THEN
        IF (XMA.LT.1.E-6.AND.XMB.LT.1.E-6) RETURN
        IF (XMA.LT.1.E-6) P=(1.-ERM)**XMB
        IF (XMB.LT.1.E-6) P=(1.-ERM)**XMA
      ELSE
        P=(((XMA+XMB)*ERM/XMA)**XMA)*(((XMA+XMB)*(1.-ERM)/XMB)**XMB)
      END IF
      IF (P.LT.RF(0)) GO TO 100
      RETURN
      END
*   RF.FOR
!======================================================================
      FUNCTION RF(IDUM)
*--generates a uniformly distributed random fraction between 0 and 1
*----IDUM will generally be 0, but negative values may be used to
*------re-initialize the seed
      SAVE MA,INEXT,INEXTP
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF/0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 50 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF (MK.LT.MZ) MK=MK+MBIG
          MJ=MA(II)
50      CONTINUE
        DO 100 K=1,4
          DO 60 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I)=MA(I)+MBIG
60        CONTINUE
100     CONTINUE
        INEXT=0
        INEXTP=31
      END IF
200   INEXT=INEXT+1
      IF (INEXT.EQ.56) INEXT=1
      INEXTP=INEXTP+1
      IF (INEXTP.EQ.56) INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ=MJ+MBIG
      MA(INEXT)=MJ
      RF=MJ*FAC
      IF (RF.GT.1.E-8.AND.RF.LT.0.99999999) RETURN
      GO TO 200
      END
!======================================================================
	subroutine cal_var

	Include 'common.txt'
	INCLUDE 'PROPERTY.TXT'

	integer:: nc, nx, i, n
	real::flow0, flowL, vels(3), vells(3), Tx, Ty, Tz, T, rho

	Rm0=0
	RmL=0
	do i=1, 2
	do nc=1, ncy
	  if(i==1) nx=1
	  if(i==2) nx=ncx
	  mc=(nc-1)*ncx+nx
	  do n = 1, 3
	     vels(n)  = cs(n+1,mc,1) / cs(1,mc,1)
	     vells(n) = cs(n+4,mc,1) / cs(1,mc,1)
	  end do
	  Tx = (sp(5,1)/Boltz)*(vells(1)-vels(1)**2)
	  Ty = (sp(5,1)/Boltz)*(vells(2)-vels(2)**2)
	  Tz = (sp(5,1)/Boltz)*(vells(3)-vels(3)**2)
	  T =(Tx+Ty+Tz)/3.0
	  den =Fnum*cs(1,mc,1)/(nsmp*cc(5))

	  p =Boltz*T*den
*	  if((i==2).and.(nc==(mcell/2))) pres=p/(pout*pref)
	  rho=den*sp(5,1)  !p_e(nc)/(gasR*T_e(nc)) 
	  !rho_e(nc)=rho_j+(pout*pref-p_j(nc))/(asound**2)		  

	  if(i==1)then
		Rm0=Rm0+rho*vels(1)*(fh/ncy)
	  else
		RmL=RmL+rho*vels(1)*(fh/ncy)
	  endif
	enddo
	enddo
	Rm0=Rm0/fh
	RmL=RmL/fh
	mid0=((ncy/2)-1)*ncx+1
	midL=((ncy/2))*ncx
	u0= cs(2,mid0,1) / cs(1,mid0,1)
	uL= cs(2,midL,1) / cs(1,midL,1)

	endsubroutine cal_var


*****************************************
*   
      SUBROUTINE WallHeatFlux
	!Subroutine for implementation of Wall heat flux for Nozzle code
	!(Adding by Akhlaghi 04-21-2013)
*
      Include 'common.txt'
	INCLUDE 'PROPERTY.TXT'
	CHARACTER (len=35)	::	fname2
*
	Open(222,file='InputData.txt')
	READ(222,*)	NSM
	Close(222)

		IF (NPR.GT.NSM) THEN
		iter=iter+1
		ENDIF
	L=1
	DTWS=1000*NIS*DTM
	WRITE (fname2,'(a,I10,a)')'WallProp-NPR=',NPR,'.plt'
	OPEN (300,file=fname2)
          WRITE (300,*)'Variables=i,CollRate,Qw,Qsm,Tw,Tsm,DTw,Etr_in,
     &	Erot_in,Etr_out,Erot_out,Vx'
	
		DO i=1,NCX
			CollRate(i)=CSS(1,i,L)/DTWS
			Et_in=CSS(5,i,L)*Fnum/(DTWS*SQRT(CG(3,i)**2+CG(6,i)**2))
			Er_in=CSS(7,i,L)*Fnum/(DTWS*SQRT(CG(3,i)**2+CG(6,i)**2))
			Et_out=CSS(6,i,L)*Fnum/(DTWS*SQRT(CG(3,i)**2+CG(6,i)**2))
			Er_out=CSS(8,i,L)*Fnum/(DTWS*SQRT(CG(3,i)**2+CG(6,i)**2))
			QQ=Et_in+Er_in+Et_out+Er_out
			
			IF (NPR.GT.NSM) THEN
				QQs(i)=QQs(i)+QQ
				TTs(i)=TTs(i)+Tp
				Qsm(i)=QQs(i)/REAL(iter)
				Tsm(i)=TTs(i)/REAL(iter)
			ElSE
				iter=0 
				DO j=1,NCX
					QQs(j)=0
					TTs(j)=0
				ENDDO 
			ENDIF


			IF ((ITypeQw==1).and.(i.GT.INZ(2))) THEN
				Tp=Tw(i)
				DTw=Tw(i)*QwRF*(QQ-DesQw(i))/ABS(10*Rmdot*Cp*Tavg/aLn)
				Tw(i)=Tw(i)+DTw
			Else
				Tp=Tw(i)
			end if 

		WRITE (300,99007)i,CollRate(i),QQ,Qsm(i),Tp,Tsm(i),DTw,Et_in,
     &			  Er_in,-Et_out,-Er_out,Vx(i)
		ENDDO
	CLOSE (300)
99007 FORMAT (I8,11ES20.4)
	
		CSS(:,:,:)=0.0
	
*
	END SUBROUTINE
*****************************************






!======================================================================
      SUBROUTINE DATA2
*
*--defines the data for a particular run of DSMC2.FOR.
*
      Include 'common.txt'
	INCLUDE 'PROPERTY.TXT'


*--set data (must be consistent with PARAMETER variables)
*	
	Open(221,file='InputData.txt')
	READ(221,*)	NSM
      READ(221,*)	
	READ(221,*)	NCX
      READ(221,*)	NCY
	READ(221,*)	
	READ(221,*)	NBX
	READ(221,*)	NBY

      READ(221,*)	NSCX
      READ(221,*)	NSCY
      READ(221,*)	IFCX
      READ(221,*)	CWRX
      READ(221,*)	IFCY
      READ(221,*)	CWRY
      READ(221,*)	IIS
      READ(221,*)	ISG
      READ(221,*)	FTMP
	READ(221,*)	PIN 
      FND=PIN/(BOLTZ*FTMP)
*--FND is the number densty
	READ(221,*)	POUT
      READ(221,*)	VFX
      READ(221,*)	VFY
      READ(221,*)	FSP(1)
      READ(221,*)	FNUM
	READ(221,*)	DTM
      READ(221,*)	CB(1)
      READ(221,*)	CB(2)
      READ(221,*)	CB(3)
      READ(221,*)	CB(4) 
	READ(221,*)	CB(5)
	CB(5)=CB(4)-CB(5)
	READ(221,*)	CB(6)
	!BUFFER INFORMATION
	CW=(CB(2)-CB(1))/NCX
	CH=(CB(4)-CB(3))/NCY
	H_BUFFER=(NBY-NCY)*CH

      CB(3)=CB(3)+H_BUFFER
      CB(4)=CB(4)+H_BUFFER	!/2. FOR SYMMETRY
	CB(5)=CB(5)+H_BUFFER	!CB(4)/1.50
	CB(6)=CB(6)+H_BUFFER		!INLET HEIGHT
	
	CB(7)=CB(2)
	CB(8)=CB(7)+NBX*CW
	READ(221,*)	CB(9)
	CB(10)=CB(4)

	READ(221,*)	INZ(1)
	READ(221,*)	INZ(2) 
	READ(221,*)	INZ(3)
	READ(221,*)	aLn

*--the simulated region is from x=XB(1) to x=XB(2)
      READ(221,*)	IB(1)
      READ(221,*)	IB(2)
      READ(221,*)	IB(3)
      READ(221,*)	IB(4)
	READ(221,*)	IB(5)
	READ(221,*)	IB(6)

      READ(221,*)	ISURF(1)
	READ(221,*)	ISURF(2)
	READ(221,*)	ISURF(3)

      READ(221,*)	LIMS(1,1)
      READ(221,*)	LIMS(1,2)
      READ(221,*)	LIMS(1,3)
      READ(221,*)	LIMS(2,1)
      READ(221,*)	LIMS(2,2)
      READ(221,*)	LIMS(2,3)
      
	READ(221,*)	TSURF(1)
      READ(221,*)	TSURF(2)
	READ(221,*)	TSURF(3)

      READ(221,*)	IJET
      READ(221,*)	SP(1,1)
      READ(221,*)	SP(2,1)
      READ(221,*)	SP(3,1)
      READ(221,*)	SP(4,1)
      READ(221,*)	SP(5,1)
	READ(221,*) Cp
      READ(221,*)	ISPR(1,1)
      READ(221,*)	SPR(1,1,1)
      READ(221,*)	ISPR(2,1)
      READ(221,*)	ISP(1)
      READ(221,*)	NIS
      READ(221,*)	NSP
      READ(221,*)	NPS
      READ(221,*)	NPT
	mfp=1.0/(SQRT(2.0)*(pi*sp(1,1)**2)*FND)
      kn=mfp/(CB(4)*2.)	 
*

	READ(221,*)	ITypeQw
	READ(221,*)	Qw1
	READ(221,*)	Tw1
	READ(221,*)	QwRF
	READ(221,*)	QwEPS
	READ(221,*)	Twb
	READ(221,*)	NTwES
	READ(221,*)	NSQ
	
	do j=1,NCX
		DesQw(j)=Qw1
		Tw(j)=Tw1
	enddo 

	If ((ITypeQw==1) .and. (NTwES==1)) then 
		Open (124,file='Tw.txt')
		Do j=1,NCX
			read(124,*) Tw(j) 
		End do 
	End IF 
	
	If (ITypeQw==2) then 
		Open (124,file='Tw.txt')
		Do j=1,NCX
			read(124,*) Tw(j) 
		End do 
	End IF 

	do j=NCX+1,NCX+10
		Tw(j)=Twb
	end do 
	Close(221)
      RETURN
      END
