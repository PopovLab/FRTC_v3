module trajectory
    use kind_module
    implicit none

    integer, parameter :: mpnt = 100000


    integer nrefj(mpnt)
    !!common/refl/nrefj(mpnt)

    integer mbeg(mpnt),mend(mpnt),mbad(mpnt)

contains

subroutine init_trajectory
    use constants
    use driver_module
    implicit none
    nrefj = 0
    
    dland = zero
    dcoll = zero
    perpn = zero 
    dalf  = zero
    vel = zero
    jrad = zero
    iww = zero
    tetai = zero
    xnpar = zero
    izz = zero

    mbeg = zero
    mend = zero
    mbad = zero

end subroutine 

subroutine view(tview, ispectr,nnz,ntet) !sav2008
!!!writing trajectories into a file
    use constants
    use approximation
    use plasma
    use dispersion_module, only: zatukh
    use rt_parameters, only :  nr, itend0, kv, nmaxm    
    use dispersion_module
    use driver_module !, only: jrad, iww, izz, length

    implicit none
    
    real(wp), intent(in) :: tview

    integer, intent(in) :: ispectr, nnz, ntet  !sav#
    !common /bcef/ ynz,ynpopq
    !common /vth/ vthc(length),poloidn(length)
    real(wp) vthcg,npoli
    !common /a0ghp/ vlf,vrt,dflf,dfrt
    
    integer i, n, itr, ntraj
    integer jrc,nturn,ib,ie,jr,ifast,idir,iv
    integer jznak,jdlt,mn,mm,jchek,itet,inz
    integer, parameter :: unit_bias = 10
    integer, parameter :: m=7
    real(wp), parameter :: pleft=1.d-10 !m may be chaged together with name(m)

    real(wp) :: htet, h, xr, xdl, xdlp, xly, xlyp, xgm, xgmp, th
    real(wp) :: x, xx, z, zz, pl, pc, pa
    real(wp) :: pdec1z, pdec3z, pintld, pintal
    real(wp) :: cotet, sitet
    real(wp) :: v, refr, dek3, parn, argum
    real(wp) :: df, powpr, powd, powal, pil, pic
    real(wp) :: powcol, pia, pt, denom, powdamped, domin, fff

    character(40) fname

    print *, 'view_time=',tview
    !print *, name(m)
    if (ispectr>0) then
        write(fname,'("lhcd/traj/pos/", f9.7,".dat")') tview
    else
        write(fname,'("lhcd/traj/neg/", f9.7,".dat")') tview
    endif
    print *, fname
    !name(m) = fname
    !print *, name(m)

    htet=zero
    h=1d0/dble(nr+1)
    if(ntet.ne.1) htet=(tet2-tet1)/(ntet-1)
    open(1,file='lhcd/out/lcms.dat')
    write(1,*)'     R(m)            Z(m)'
    write(1,*)
    xr=1.d0
    xdl=fdf(xr,cdl,ncoef,xdlp)
    xly=fdf(xr,cly,ncoef,xlyp)
    xgm=fdf(xr,cgm,ncoef,xgmp)
    do i=1,101
        th=dble(i-1)*pi2/dble(100)
        cotet=dcos(th)
        sitet=dsin(th)
        xx=-xdl+xr*cotet-xgm*sitet**2
        zz=xr*xly*sitet
        x=(r0+rm*xx)/1d2
        z=(z0+rm*zz)/1d2
        write(1,5) x,z
    end do
    close(1)

    open(1,file=fname)
    write(1,3) !write header 
    ntraj=0 !sav2008
    do itr=1,nnz*ntet !sav2008
        pow=1.d0
        pl=zero
        pc=zero
        pa=zero
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
        jrc=nr+1
        jznak=-1
        nturn=1
        if(mbad(itr).eq.0) then 
            ntraj=ntraj+1
            ib=mbeg(itr)
            ie=mend(itr)
10          continue
            do i=ib,ie
                v=vel(i)
                jr=jrad(i)
                refr=perpn(i)
                npoli=poloidn(i)
                ifast=iww(i)
                vthcg=vthc(i)
                idir=izz(i)
                dek3=zero
                th=tetai(i)
                parn=xnpar(i)
                if(itend0.gt.0) then
                    argum=clt/(refr*valfa)
                    dek3=zatukh(argum,abs(jr),vperp,kv)
                end if

                call distr(v,abs(jr),iv,df)
                if(jr.lt.0) then    !case of turn
                    jr=-jr
                    !variant          pintld=-dland(i)*df
                    !!          pintld=-dland(i)*(dflf+dfrt)/2d0
                    pintld=dabs(dland(i)*(dflf+dfrt)/2d0)
                    pdec2=dexp(-2d0*dcoll(i))
                    pintal=dabs(dalf(i)*dek3)
                else
                    pdec2=dcoll(i)
                    pdecv=dland(i)
                    !!          pdec1=-pdecv*df
                    pdec1=dabs(pdecv*df)
                    pdec3=dabs(dalf(i)*dek3)
                    pintld=(pdec1+pdec1z)/2d0*h
                    pintal=(pdec3+pdec3z)/2d0*h
                    pdec1z=pdec1
                    pdec3z=pdec3
                end if
                powpr=pow
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
                pil=pintld
                pic=.5d0*dabs(dlog(pdec2))
                pia=pintal
                pt=1.d0-pow  !total absorbed power
                denom=pil+pic+pia
                powdamped=1.d0-dexp(-2.d0*denom)
                domin=powpr*powdamped
                if(denom.ne.zero) then
                    fff=domin/denom
                    pl=pl+dabs(pil*fff)  !el. Landau absorbed power
                    pc=pc+dabs(pic*fff)  !el. collisions absorbed power
                    pa=pa+dabs(pia*fff)  !alpha Landau absorbed power
                end if
                xr= rho(i)!h*dble(jr)
                cotet=dcos(th)
                sitet=dsin(th)
                xdl=fdf(xr,cdl,ncoef,xdlp)
                xly=fdf(xr,cly,ncoef,xlyp)
                xgm=fdf(xr,cgm,ncoef,xgmp)
                xx=-xdl+xr*cotet-xgm*sitet**2
                zz=xr*xly*sitet
                x=(r0+rm*xx)/1d2
                z=(z0+rm*zz)/1d2 
                jdlt=jr-jrc
                jrc=jr
                if(jdlt*jznak.lt.0.and.nturn.lt.m-1) then
                    nturn=nturn+1
                    jznak=-jznak
                end if

                        !   R, Z, rho, theta, N_par, N_pol, P_tot,P_land, P_coll, vth,  'slow=1','out=1', N_traj 
                write(1, 7) x, z, xr,  th,    parn,  npoli, pt,   pl,     pc,     vthcg, ifast,  idir,    itr
                
                if(pt.ge.1d0-pleft) go to 11 !maximal absorbed power along a ray
            end do
            jchek=jrad(ie+1)
            if(jchek.ne.0) then  !continue this trajectory
                ib=idnint(dland(ie+1))
                ie=idnint(dcoll(ie+1))
                goto 10
            end if
11          continue

        end if
    end do


    close(1)


1     format(2x,'N_traj',3x,'mbad',6x,'theta',9x,'Npar',9x,'rho_start')
2     format('R_pass',4x,'Ptot',6x,'Pland',6x,'Pcoll',8x,'Pa',7x,'dPtot',6x,'dPland',5x,'dPcoll',6x,'dPa')
3     format(5x,'R',10x,'Z',11x,'rho',8x,'theta',7x,'N_par',7x,'N_pol',6x,'P_tot',7x,'P_land',6x,'P_coll',6x,'vth',4x,'slow=1',4x,'out=1',2x,'N_traj',6x)
4     format(i3,5x,8(f6.3,5x))
5     format(6(e13.6,3x))
6     format(2(i6,2x),4(e13.6,1x))
7     format(10(e11.4,1x),i5,2x,i5,2x,i5)
8     format('after radial pass=',i3,2x,' P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
9     format('Total passes:           P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
20    format('written time slice (seconds) =',f9.3)
    end

    subroutine traj(xm0, tet0, xbeg, nmax, nb1, nb2, nomth, nomnz, pabs) !sav2009
        use constants, only : tiny1
        use rt_parameters, only: eps, rrange, hdrob, nr, ipri, iw
        use dispersion_module, only: iroot, ider, izn
        use dispersion_module, only: xnr1, xnr2, xnr3, xnr4
        use dispersion_module, only: extd4, disp2, disp2_iroot3
        use driver_module, only: inak, im4, hrad, irs, iabsorp, iznzz, iwzz, irszz, rzz
        use driver_module, only: tetzz, xmzz
        use driver_module, only: driver2, driver4
        use dispersion_equation, only: ynz
        implicit none
        real(wp), intent(in)    :: xm0
        real(wp), intent(in)    :: tet0
        real(wp), intent(inout) :: xbeg
        real(wp), intent(in)    :: pabs
        integer,  intent(inout) :: nmax        
        integer,  intent(inout) :: nb1, nb2        
        integer,  intent(in)    :: nomth, nomnz

        integer :: nrefl
        integer :: irep
        integer :: irf, irf1
        integer :: ib2
        integer :: irs0        
        integer :: inak_saved 
        real(wp), parameter :: pgdop=0.02d0
        real(wp), parameter :: hmin=0.d-7 !sav2008, old hmin=1.d-7
        real(wp) :: eps0
        real(wp) :: rrange0, hdrob0, tet, xm, hr
        real(wp) :: xsav, xend,hsav, h1
        real(wp) :: ystart(2),yy(4)

        real(wp) :: xnr, prt, prm
        real(wp) :: ynz0, x1, x2, rexi, tetnew
        real(wp) :: xmnew, rnew, xnrnew, xnrv
        real(wp) :: pg1, pg2, pg3, pg4, pg

        print *, '----------------- traj ----------------'
        eps0=eps
        rrange0=rrange
        hdrob0=hdrob
        iroot=1
        nrefl=0
        ider=1
        im4=0
        nb1=0
        nb2=0
        irep=0
        tet=tet0
        xm=xm0
        hr=1.d0/dble(nr+1) !sav2008
        hrad=hr

        !---------------------------------------
        ! find saving point and define
        ! parameters, depending on direction
        !---------------------------------------
  
  10    print *, '---------- label 10 ----------'
        print *, '  hr=', hr, 'irs=', irs
        print *, 'xbeg=', xbeg
        irf1=idnint(xbeg/hr)
        if (dabs(irf1*hr-xbeg).lt.tiny1)  then
            xsav=hr*irf1
            print *, 'xxxx'
        else
            print *, 'yyyy irs=', irs
            irf=int(xbeg/hr)
            print *, 'xsav= ', xsav
            print *, 'irf= ', irf
            if (irs.eq.1)  xsav=hr*irf
            if (irs.eq.-1) xsav=hr*(irf+1)
        end if
        xend= 0.5d0 - 0.5d0*irs + tiny1*irs
        if (ipri.gt.2) write (*,*) 'xbeg-xend',xbeg,xend
        hsav = -hr*irs
        h1 = hsav
        !---------------------------------------
        ! solve eqs. starting from xbeg
        !---------------------------------------
        ystart(1) = tet
        ystart(2) = xm
        
        print *, '   h1=', h1
        print *, ' xbeg=', xbeg
        print *, ' xend=', xend
        print *, ' xsav=', xsav
        print *, '------'
        call driver2(ystart,xbeg,xend,xsav,hmin,h1, pabs)
        tet = ystart(1)
        xm = ystart(2)
        ib2 = 0

        !---------------------------------------
        ! absorption
        !---------------------------------------
        if(iabsorp.ne.0) then
            if(ipri.gt.2) write (*,*)'in traj() iabsorp=',iabsorp
            nmax=nrefl
            return
        end if
        if (xend.eq.xbeg) nb1=nb1+1
        !sav2008 20    continue
  
        !--------------------------------------------------------
        !  pass turning point
        !--------------------------------------------------------
        irs0=irs
        ider=0
        print *, '--- pass turning point ---'
        call disp2(xend,xm,tet,xnr,prt,prm)
        ider=1
        ynz0 = ynz
  40    yy(1)=tet 
        yy(2)=xm
        yy(3)=xend
        yy(4)=xnr
        x1=0d0
        x2=1d+10
        rexi=xend
        inak_saved = inak
        print *, '---- driver4 -----'
        call driver4(yy,x1,x2,rexi,hmin, extd4)
        if(iabsorp.eq.-1) return !failed to turn

        tetnew = yy(1)
        xmnew  = yy(2)
        rnew   = yy(3)
        xnrnew = yy(4)

        if(ipri.gt.2) write (*,*) 'from r=',rexi,'to r=',rnew
        print *, 'from r=',rexi,'to r=',rnew
        print *, '   tet =', tet
        print *, 'tetnew =', tetnew
        print *, ' dteta =',  tetnew - tet
        if (abs(tetnew - tet)> 1.0) then
            pause 'bug dt'
        endif
        !---------------------------------------
        ! find mode
        !---------------------------------------
        iroot=3
        ider=0
        xnrv=xnrnew
        call disp2_iroot3(rnew,xmnew,tetnew,xnrv,prt,prm)
        ider=1
        iroot=1
        !ipric      if (ipri.gt.2) then
        !ipric       write (*,*)'nr check, r=',rnew,' tet=',tetnew
        !ipric       write (*,*)'iw=',iw,' izn=',izn
        !ipric       write (*,*) xnrnew,xnr1
        !ipric       write (*,*) xnr2,xnr3,xnr4
        !ipric       pause
        !ipric      end if
        pg1 = abs(xnrnew-xnr1)
        pg2 = abs(xnrnew-xnr2)
        pg3 = abs(xnrnew-xnr3)
        pg4 = abs(xnrnew-xnr4)

        pg = dmin1(pg1,pg2,pg3,pg4)
        if (dabs(pg/xnrnew).gt.pgdop) then
            !---------------------------------------------
            ! bad accuracy, continue with 4 equations
            !--------------------------------------------
            ib2=ib2+1
            nb2=nb2+1
            if (ib2.gt.4) then
                if (ipri.gt.1) write (*,*) 'error: cant leave 4 eqs'
                iabsorp=-1
                return
            end if
            eps=eps/5d0
            rrange=rrange*2d0
            hdrob=hdrob*2d0
            inak = inak_saved ! костыль - восстанавливаю занчение счетчика точек траектории
            goto 40
        end if
        !-------------------------------------
        !          change wave type
        !-------------------------------------
        if (pg.ne.pg1) then
            if (pg.eq.pg2) izn=-izn
            if (pg.eq.pg3) iw=-iw
            if (pg.eq.pg4) iw=-iw
            if (pg.eq.pg4) izn=-izn
        end if
        if (irs0.ne.irs) nrefl=nrefl+1
        xbeg=rnew
        tet=tetnew
        xm=xmnew
        im4=1
        eps=eps0
        rrange=rrange0
        hdrob=hdrob0
        print *, '---- goto 10 '
        if(nrefl.lt.nmax) goto 10
        rzz=xbeg
        tetzz=tet
        xmzz=xm
        iznzz=izn
        iwzz=iw
        irszz=irs
    end  


    

end module trajectory