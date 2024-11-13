      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      include 'SMAAspUserSubroutines.hdr'
C
      DIMENSION TIME(2)
C
      IF(LOP == 0) THEN
C Initialization in UEXTERNALDB/VEXTERNALDB
          call MutexInit(1)
          call MutexInit(2)
          ! initialize Mutex #1
      end if

      RETURN
      END


      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS(3)
      CHARACTER*80 SNAME

      real(kind=8) x, y, z
      real(kind=8) P(3), startP(5), endP(5), v(3)
      real(kind=8) delta(3)

      ! 声明一个动态二维数组, 用于存储文件中的数据
      real(kind=8), dimension(:,:), allocatable :: P_array
      ! 声明一些变量
      integer :: i, j, nrow, ncol
      character (len=500) line ! 每行的内容
      character (len=500) fn
      data iread /1/
      save iread
      save P_array,nrow,ncol

      if (iread == 1) then
          iread = 2
          fn='D:\SIMULIA\temp\tubularStructure\Path.txt'
C     读取路径文件
          ! 打开文件
          open(unit=19, file=fn, status='old', action='read', iostat=i)
          if (i /= 0) then
              print *, 'Error opening file'
              stop
          end if

          ! 初始化行数和列数
          nrow = 0
          ncol = 5

          ! 循环读取每一行
          do
              read(19, '(A)', end=100) line
              nrow = nrow + 1 ! 行数加一
          end do
100       continue
          rewind(19)

          ! 分配动态数组的空间
          allocate(P_array(nrow, ncol))

          ! 从文件中读取数据并赋值给动态数组
          do i = 1, nrow
              read (19, *) (P_array(i,j), j = 1, ncol)
          end do
          ! 关闭文件
          close(19)
          print*,"Path file has read successfully"
          print*,iread
      end if

C     1. 判断当前积分时刻所处的时间段。
      t = TIME(2)
      do i = 1,nrow-1
          startP = P_array(i,:)
          endP = P_array(i+1,:)
          if(t <= endP(1))then
              exit
          end if
      end do
C     2. 获取该段路径的起始时间，并计算该时间段内的速度
      t0 = startP(1)
      ! 确定速度大小及方向
      do i=1,3
          v(i) = (endP(i+1) - startP(i+1))/(endP(1) - startP(1))
      end do
C     3. 计算当前时刻的加载位置
      do i = 1,3
          if (v(i) /= 0.0) then
              P(i) = v(i)*(t - t0) + startP(i+1) - v(i)/abs(v(i))*0.3
              delta(i) = (0.35 + 0.45) / 2.0
          else
              P(i) = v(i)*(t - t0) + startP(i+1)
              delta(i) = 0.3
          end if
      end do

      x = COORDS(1)
      y = COORDS(2)
      z = COORDS(3)

C 确定载荷施加区域    
C startP(5) 表示载荷状态，当startP(5)=0时表示空行程，没有压力
      if((abs(x - P(1)) < delta(1)) .and. (abs(y - P(2)) < delta(2))
     1 .and. (abs(z - (P(3) - 0.2)) < 0.025)) then
          F = 0.152 * startP(5)
      else
          F = 0.0
      endif

      RETURN
      END


      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      include 'SMAAspUserSubroutines.hdr'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION DSTRES(6),D(3,3)

      REAL(kind=8) Lambda(16)
      REAL(kind=8) E1(16),E2(16),E3(16),G12(16),G13(16),G23(16)
      REAL(kind=8)::h(16)=0.0
      REAL(kind=8)::R(16) = 0.0
      real(kind=8) B(5)                       !用于存储移位函数(Shift function)的拟合参数
      real(kind=8)::cTM = 0.0                 !用于存储C_ij所对应的移位函数cTM_ij
      real(kind=8)::Eng(6,6,16) = 0.0         !用于存储工程常数
      real(kind=8)::f(6,6) = 0.0              !用于存储刚度系数
      real(kind=8)::C0(6,6) = 0.0             !用于存储刚度矩阵
      real(kind=8)::C(6,6,16) = 0.0           !用于存储刚度矩阵
      real(kind=8)::s(6,16) = 0.0
      real(kind=8) E1_inf,E2_inf,E3_inf,G12_inf,G23_inf,G13_inf
      real(kind=8) v12,v13,v23,tau
      real(kind=8) E1w, E2w, E3w
      real(kind=8) TM, Px, t1
      INTEGER tempvar
      INTEGER w
      logical :: firstrun = .true.
C debug
      !IF (firstrun) then
      !    write (*,*)"Please input an integer:"
      !    read(*,*)tempvar
      !    firstrun = .false.
      !end if
      !tempvar=1234

C Initialization variables
C N: The numbers of MaxWell element
C B1 B2 B3 —— Parameters of Shift function equation
      call MutexLock(1)
      N = 16
      B = (/6.84465752463272,-0.147327023748936,0.000342094310052091,
     1 6.38330404062246E-6,-3.17644160305058E-8/)
	 
C Calculate the median temperature and WLF coefficient
      TM = TEMP !+ DTEMP/2.0
      Px = B(1) + B(2)*TM + B(3)*TM**(2.0) + B(4)*TM**(3.0)
     1 + B(5)*TM**(4.0)
      cTM = 10.0**(-Px)
	  
C Define relaxation time and modulus.
      Lambda = (/0.0001,0.001,0.01,0.1,1.0,10.0,100.0,1000.0,10000.0,
     1 100000.0,1000000.0,10000000.0,100000000.0,1000000000.0,
     2 10000000000.0,100000000000.0/)
      E1 = (/0.00852,208.70905,445.91378,243.8356,181.40715,95.22916,
     1 49.70814,96.48326,119.03913,79.68588,88.59901,44.68303,73.09294,
     2 20.34076,44.12642,23.32033/)


      E2 = (/99.875,350.69214,132.75968,160.84225,91.75791,58.79377,
     1 36.16295,51.51099,53.51274,52.52625,24.81714,33.62895,20.45954,
     2 36.37052,7.315,41.80219/)

      E3 = (/109.70815,225.09451,137.8043,100.42143,75.12228,44.3639,
     1 34.20719,46.0066,30.16637,29.49191,18.00616,29.82038,24.15735,
     2 21.69017,2.73976,7.66513/)

      G12 = (/68.76763,141.09449,86.37895,62.94649,47.08839,27.80832,
     1 21.44186,28.83801,18.90899,18.48622,11.28668,18.69211,15.14239,
     2 13.59591,1.71735,4.80468/)

      G13 = (/44.42965,91.15886,55.80804,40.66871,30.42305,17.96651,
     1 13.85324,18.63177,12.21679,11.94365,7.29214,12.07667,9.78325,
     2 8.7841,1.10955,3.10423/)

      G23 = (/47.47431,97.40578,59.63244,43.45565,32.50787,19.19771,
     1 14.80257,19.90857,13.05398,12.76212,7.79186,12.90425,10.45368,
     2 9.38605,1.18558,3.31695/)
	  
      E1_inf = 30.0
      E2_inf = 20.0
      E3_inf = 10.0
      G12_inf = 6.26823
      G13_inf = 4.0498
      G23_inf = 4.32733
      
      IF (TM.GE.220.0) THEN
          DO w = 1, N
              E1(w) = E1(w)/1844.18 * 30.0
              E2(w) = E2(w)/1272.83 * 20.0
              E3(w) = E3(w)/946.47 * 10.0 
              G12(w) = G12(w)/593.27 * 6.27
              G13(w) = G13(w)/383.3 * 4.05
              G23(w) = G23(w)/409.57 * 4.33
          END DO
          E1_inf = E1_inf/1844.18 * 30.0
          E2_inf = E2_inf/1272.83 * 20.0
          E3_inf = E3_inf/946.47 * 10.0 
          G12_inf = G12_inf/593.27 * 6.27
          G13_inf = G13_inf/383.3 * 4.05
          G23_inf = G23_inf/409.57 * 4.33
      END IF
	  
C 将工程常数写入Eng矩阵
      Eng(1,1,:) = E1
      Eng(1,2,:) = E1
      Eng(1,3,:) = E1
      Eng(2,1,:) = E2
      Eng(2,2,:) = E2
      Eng(2,3,:) = E2
      Eng(3,1,:) = E3
      Eng(3,2,:) = E3
      Eng(3,3,:) = E3
      Eng(4,4,:) = G12
      Eng(5,5,:) = G13
      Eng(6,6,:) = G23


      t1 = TIME(2)*cTM
      call calc_E(t1,E1,Lambda,N,E1_inf,E1w)
      call calc_E(t1,E2,Lambda,N,E2_inf,E2w)
      call calc_E(t1,E3,Lambda,N,E3_inf,E3w)
      
      v12 = 0.41
      v13 = 0.32
      v23 = 0.40
      v21 = v12*E2w/E1w
      v31 = v13*E3w/E1w
      v32 = v23*E3w/E2w
      tau = 1.0- v12*v21 - v23*v32 - v31*v13 - v21*v32*v13 - v21*v23*v31
C 计算刚度系数矩阵      
      f(1,1) = (1.0 - v23*v32)/tau
      f(1,2) = (v21 + v31*v23)/tau
      f(1,3) = (v31 + v21*v32)/tau
      f(2,1) = (v12 + v32*v13)/tau
      f(2,2) = (1.0 - v13*v31)/tau
      f(2,3) = (v32 + v12*v31)/tau
      f(3,1) = (v13 + v12*v23)/tau
      f(3,2) = (v23 + v21*v13)/tau
      f(3,3) = (1.0 - v12*v21)/tau
      f(4,4) = 1.0
      f(5,5) = 1.0
      f(6,6) = 1.0

C Calculate characteristic function
      do w = 1,N
          R(w) = EXP(-(cTM*DTIME)/Lambda(w))
      end do

      DO w = 1, N
          h(w)=(1.0/DTIME)*(Lambda(w)/cTM)*(1.0-R(w))
      END DO
      
C 计算刚度矩阵
      DDSDDE = 0.0
      C0(1,1) = f(1,1)*E1_inf
      C0(1,2) = f(1,2)*E1_inf
      C0(1,3) = f(1,3)*E1_inf
      C0(2,1) = f(2,1)*E2_inf
      C0(2,2) = f(2,2)*E2_inf
      C0(2,3) = f(2,3)*E2_inf
      C0(3,1) = f(3,1)*E3_inf
      C0(3,2) = f(3,2)*E3_inf
      C0(3,3) = f(3,3)*E3_inf
      C0(4,4) = f(4,4)*G12_inf
      C0(5,5) = f(5,5)*G13_inf
      C0(6,6) = f(6,6)*G23_inf

C 初始化刚度矩阵      
      do i = 1,6
          do j = 1,6
              do w = 1,N
                  C(i,j,w) = f(i,j)*Eng(i,j,w)
              end do
              DDSDDE(i,j) = C0(i,j)
          end do
      end do

C 计算刚度矩阵      
      do i = 1,6
          do j = 1,6
              do w = 1,N
                  DDSDDE(i,j) = DDSDDE(i,j) + C(i,j,w)*h(w)
              end do
          end do
      end	do
      
      s = reshape(STATEV(7:102), [6,16])
C Relaxed and unrelaxed stress at initialization t=0
      IF (TIME(2).EQ.0.0) THEN
          DO i = 1,NTENS
              STATEV(i) = 0.0
              !STRESS_rel(i)
              DO w = 1, N
                  s(i,w) = 0.0
                  !STRESS_unrel(i,w)
              END DO
          END DO
      END IF

C Calculate Stress 
C     σ_i=σ_i0  + DDCDDE*DSTRAN+Rw*σ_iw(n-1)  
      DO i = 1, NTENS
          STRESS(i) = STATEV(i)
          !STRESS(i)=STRESS(i0)
          DO j = 1, NTENS
              STRESS(i) =  STRESS(i) + DDSDDE(i,j)*DSTRAN(j)
              !STRESS(i)=STRESS(i0)+DDSDDE(i,j)*DSTRAN(j)
          END DO
          DO w = 1, N
              STRESS(i) = STRESS(i) + R(w)*s(i,w)
          END DO
      END DO

C 更新应力分量  σ_i0      
      DO i = 1,NTENS
          DO j = 1, NTENS
              STATEV(i) = STATEV(i) + C0(i,j)*DSTRAN(j)
          END DO
      END DO

      DO i = 1,NTENS
          do w = 1, N
C Calculate the initial value of STRESS_unrel when t=t_n  
              s(i,w) = R(w)*s(i,w)
C Add the contribute of Cijw
              DO j = 1,NTENS
                  s(i,w) = s(i,w) + C(i,j,w)*DSTRAN(j)*h(w)
              ENDDO
          END DO
      END DO
      STATEV(7:102) = reshape(s,[96])
      call MutexUnlock(1)
      RETURN
      END


      SUBROUTINE calc_E(t, Ew, lambda, N, Einf, E)
      include 'SMAAspUserSubroutines.hdr'
      REAL(kind=8), INTENT(IN) :: t, Ew(N), lambda(N), Einf
      REAL(kind=8), INTENT(OUT) :: E

      INTEGER :: w
      call MutexLock(2)

      E = Einf

      DO w = 1, N
          E = E + Ew(w)*EXP(-t/lambda(w))
      END DO

      call MutexUnlock(2)
      END SUBROUTINE calc_E

