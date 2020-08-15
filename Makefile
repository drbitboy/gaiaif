include Makefile.dotar
	

Makefile.dotar:
	echo "tar:@#tar zcf - 00readme.txt *.m *.py Makefile | tee gaiaif.tar.gz | tar zdvf -" | tr '@#' \\n\\t > $@

test:
	[ -d "jsonlab/" ] || git clone https://github.com/fangq/jsonlab.git
	echo 'addpath("jsonlab/") ; x.fov={{1,2},3} ; x.limit=2 ; x.j2000=1 ; x.ppm=0 ; x.mags=1 ; x.heavy=1 ; fov_cmd(x)' | octave --no-gui
