include Makefile.dotar

OBSPOS=[0 0 1e9]
OBSVEL=[7 8 9]
OBSY=2020-07-02T12:00:00.000
LIMIT=2
FOV={{1,2},2.9}


Makefile.dotar:
	echo "tar:@#tar zcf - 00readme.txt *.m *.py Makefile | tee gaiaif.tar.gz | tar zdvf -" | tr '@#' \\n\\t > $@

test:
	[ -d "jsonlab/" ] || git clone https://github.com/fangq/jsonlab.git
	( echo 'addpath("jsonlab/") ;' \
	 ; echo 'x.fov=$(FOV) ;' \
	 ; echo 'x.limit=$(LIMIT) ;' \
	 ; echo 'x.j2000=1 ;' \
	 ; echo 'x.ppm=0 ;' \
	 ; echo 'x.mags=1 ;' \
	 ; echo 'x.heavy=1 ;' \
	 ; echo 'x.obspos=$(OBSPOS) ;' \
	 ; echo 'x.obsvel=$(OBSVEL) ;' \
	 ; echo 'x.obsy="$(OBSY)" ;' \
	 ; echo 'format longG ;' \
	 ; echo 'fov_cmd(x)' \
	 ) \
	 | octave --no-gui
