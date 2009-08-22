PROGRAM correlograms
	IMPLICIT NONE
	
	! correlograms.f90: Written by Evan Thilo, University of Washington Department of Applied Mathematics, 8/21/2009.
	!				    Edited by Nicholas Cain, University of Washington Department of Applied Mathematics, 8/21/2009.
	! Please direct questions/comments to ethilo@gmail.com
	!									  nicain.seattle@gmail.com
	
	character(len=220) hstname
	character(len=220) dataFile1
	character(len=220) dataFile2

	integer(kind=4) maxspikes
	integer(kind=4) histogram
    parameter(histogram=6001)
    parameter(maxspikes=50000)
      
    integer(kind=4) auto1(histogram),auto2(histogram)
    integer(kind=4) xcor1(histogram),xcor2(histogram)
    integer(kind=4) ref,tpadstartrank1, tpadstartrank2
    integer(kind=4) i,k,d,ii
    integer(kind=4) Aprobe,Xprobe
    
    double precision corrwidth
	integer numberOfSpikes1, numberOfSpikes2
	double precision numberOfSpikes1p, numberOfSpikes2p
	double precision osc1(maxspikes)
	double precision osc2(maxspikes)

996 format(4I10)
  
    !****************************************************************************
    !**************************Auto/Cross Correlations***************************
    !****************************************************************************
      
    ! Name the output file:
    !	It will include 4 columns of bins, auto and cross correlograms,
    !	using each of the neurons as references to the other.
    !	(basically 2 sets of auto/crosses that are essentially 'inverses' of each other.)
    hstname='correlogramsOut.dat'
	
	! Name of input files:
	!	Should be a .txt file with a single column of spike times.
	dataFile1='crossStudy1_a.dat'
	dataFile2='crossStudy1_b.dat'
	
	! Read in the input data:
	open(unit = 11, file = dataFile1)
		read (11,*) numberOfSpikes1p
		numberOfSpikes1=int(numberOfSpikes1p)

		do ii=1,numberOfSpikes1
			read(11, *) osc1(ii)
		end do
	close(11)
	
	open(unit = 11, file = dataFile2)
		read (11,*) numberOfSpikes2p
		numberOfSpikes2=int(numberOfSpikes2p)

		do ii=1,numberOfSpikes2
			read(11, *) osc2(ii)
		end do
	close(11)	

	numberOfSpikes2=int(numberOfSpikes2)
	
    tpadstartrank1=1
    tpadstartrank2=1
      
    !this defines the largest ISI to include.
    corrwidth=(histogram-1)/2
      
    open(unit=13,file=hstname,status='replace')
      
	do 126 i = 1,size(auto1)
		auto1(i)=0
		xcor1(i)=0
		auto2(i)=0
		xcor2(i)=0
126  end do

	Aprobe=0
	Xprobe=0

	do 103 ref = tpadstartrank1,numberOfSpikes1
	!Autocor of osc1
	!Set the starting corr position
		
		do while (abs(osc1(ref)-osc1(Aprobe)) .gt. corrwidth)
			Aprobe=Aprobe+1
		end do
	!walk down spiketimes, record the ISI

		do 104 k = Aprobe,numberOfSpikes1
			if (abs(osc1(ref)-osc1(k)) .lt. corrwidth) then
				d=dnint(osc1(ref)-osc1(k))+corrwidth
				auto1(d)=auto1(d)+1
			else
				exit
			end if
104     end do
			
		!Xcor, osc1 as ref
		do while (abs(osc1(ref)-osc2(Xprobe)) .gt. corrwidth)
			Xprobe=Xprobe+1
		end do

		do 105 k = Xprobe, numberOfSpikes2
			if (abs(osc1(ref)-osc2(k)) .lt. corrwidth) then
				d=dnint(osc1(ref)-osc2(k))+corrwidth
				xcor1(d)=xcor1(d)+1
			else
				exit
			end if
105     end do

103  end do

	Aprobe=0
	Xprobe=0
	do 106 ref = tpadstartrank2,numberOfSpikes2
		!Autocor of osc2
		!Set the starting corr position
		do while (abs(osc2(ref)-osc2(Aprobe)) .gt. corrwidth)
			Aprobe=Aprobe+1
		end do

		!walk down spiketimes, record the ISI
		do 107 k = Aprobe,numberOfSpikes2
			if (abs(osc2(ref)-osc2(k)) .lt. corrwidth) then
				d=dnint(osc2(ref)-osc2(k))+corrwidth
				auto2(d)=auto2(d)+1
			else
				exit
			end if
107      end do
		!Xcor, osc2 as ref
		!Set the starting corr position
		do while (abs(osc2(ref)-osc1(Xprobe)) .gt. corrwidth)
			Xprobe=Xprobe+1
		end do

		!walk down spiketimes, record the ISI
		do 108 k = Xprobe, numberOfSpikes1
			if (abs(osc2(ref)-osc1(k)) .lt. corrwidth) then
				d=dnint(osc2(ref)-osc1(k))+corrwidth
				xcor2(d)=xcor2(d)+1
			else
				exit
			end if
108      end do
106  end do

	write(13,*) 'auto1, xcor1, auto2, xcor2'

	!Record all 4 correlograms to the hstgram file
	do 109 i=1,size(auto1)
		write(13,996) auto1(i),xcor1(i),auto2(i),xcor2(i)
109  end do

END PROGRAM correlograms
      