include Makefile.dotar
	

Makefile.dotar:
	echo "tar:@#tar zcf - 00readme.txt *.m *.py Makefile | tee gaiaif.tar.gz | tar zdvf -" | tr '@#' \\n\\t > $@

test:
	[ -d "jsonlab/" ] || git clone https://github.com/fangq/jsonlab.git
	@( echo 'addpath("jsonlab/") ;' \
	 ; echo 'x.fov={{1,2},2.9} ;' \
	 ; echo 'x.limit=2 ;' \
	 ; echo 'x.j2000=1 ;' \
	 ; echo 'x.ppm=0 ;' \
	 ; echo 'x.mags=1 ;' \
	 ; echo 'x.heavy=1 ;' \
	 ; echo 'x.obspos=[4 5 6] ;' \
	 ; echo 'x.obsvel=[7 8 9] ;' \
	 ; echo 'x.obsy="2020-08-17T12:34:56.789" ;' \
	 ; echo 'fov_cmd(x)' \
	 ) \
	 | octave --no-gui
