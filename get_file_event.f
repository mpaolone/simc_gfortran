      subroutine get_file_event(th_spec_e,th_spec_p,
     >    dxdz,dydz,e_mom,dxdzp,dydzp,p_mom)
c
c  input variables:
c        th_spec_e : central spec angle for electron (rad)
c        th_spec_p : central spec angle for proton (rad)
c  output variables:
c        dxdz : xptar for electron
c        dydz : yptar for electron
c        e_mom : electron momentum ( MeV)
c        dxdzp : xptar for proton
c        dydzp : yptar for proton
c        p_mom : proton momentum ( MeV)
c
      implicit none
c
         real*8 e_mom
         real*8 p_mom
         real*8 th_spec_e
         real*8 th_spec_p
         real*8 dxdz,dydz
         real*8 dxdzp,dydzp
         character*80 multpifile
         integer count,count_miss
         logical first
         logical end_of_2pi_file
         data first /.true./
         common /eventfile/  end_of_2pi_file
c
        if ( first) then
           first = .false.
           count = 0.
           count_miss = 0.
           write(*,*) ' Which input file ?'
           read(*,'(a80)') multpifile
           write(*,*) ' Opening file ',multpifile
           open(unit=51,file=multpifile)
        endif
c
c
         end_of_2pi_file = .false. 
         read(51,*,end=999,err=999) dxdz,dydz,e_mom,dxdzp,dydzp,p_mom
         count = count + 1
c
         p_mom = p_mom *1000.
         e_mom = e_mom *1000.
c        
         return
c
 999     write(*,*) ' reached end of file'
         write(*,*) ' count = ',count
         write(*,*) ' count_miss = ',count_miss,float(count-count_miss)/float(count)
         end_of_2pi_file = .true.
         return
         end
