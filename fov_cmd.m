%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fov_cmd.m
%%%
%%% Usage (see also Sample Usage below):
%%%
%%%   -- Required --
%%%
%%%   arg.fov={{1,2},{2,3},{2,1},...};  % Polygonal FOV with RA,Dec vertices, deg
%%%   arg.fov={{1,2},1};                % Circular FOV, center at RA=1deg, Dec=2deg; Radius =1deg;
%%%   arg.fov={{1,2},{2,3}};            % RA,Dec box (not a real FOV)
%%%                                     % N.B. {x,y,z} may also be used instead of any {RA,Dec}
%%%                                     % N.B. ICRS (default) or J2000l see arg.j2000 below
%%%
%%%   -- Optional --
%%%
%%%   arg.limit=200;                          % Limit to 200 stars
%%%   arg.magtype='g';                        % Use this magnitude type; g (default) or bp or rp
%%%   arg.magmax=15;                          % Limit Magnitudes to max; see .magtype
%%%   arg.magmin=11;                          % Limit Magnitudes to min; see .magtype
%%%   arg.j2000=<anything>                    % Vertices are in J2000, not ICRS (not yet implemented)
%%%   arg.buffer=0.001                        % Add buffer around FOV,degrees (not yet implemented)
%%%   arg.ppm=1                               % Include Parallax and Proper Motions in returned data
%%%   arg.mags=1                              % Include all phot_*_mean_mag values in returned data
%%%   arg.heavy=1                             % Include errors and corr. coeffs, *_error, *_corr
%%%   arg.obspos=[X Y Z]                      % Observer position vector, km (Note 1)
%%%   arg.obsvel=[VX VY VZ]                     % Observer velocity vector, km (Note 2)
%%%   arg.obsy=YYYY.ddddd                     % Observer time, fractional year (Note 3)
%%%
%%%   arg.gaiasqlite3='dirpath/gaia.sqlite3'  % Use this Gaia SQLite3 file
%%%                                           % N.B. Must end with .sqlite3
%%%                                           % Also *_heavy.sqlite3 if arg.heavy is present
%%%
%%% Notes
%%%
%%% 1) arg.obspos is solar system barycentric-relative observer
%%%    position; its presence triggers parallax correction
%%%
%%% 2) arg.obsvel is solar system barycentric-relative observer
%%%    velocity; its presence triggers stellar aberration
%%%    correction
%%%
%%% 3) arg.obsy is observer time, as a fractional year; it presence
%%%    triggers Proper Motion correction e.g. 2015.5 => 2020-06-01
%%%
%%% Returns
%%%
%%%   - array of per-star data as objects, sorted by increasing mean magnitude
%%%
%%%     - arr{index}.mean_mag   % Mean magnitude of type requested (default = g)
%%%     - arr{index}.ra         % Nominal RA of star, deg, 2015.5 epoch
%%%     - arr{index}.dec        % Nominal Dec of star, deg, 2015.5 epoch
%%%     - arr{index}.offset     % Database offset, integer
%%%     - arr{index}.parallax   % Parallax, mas, if requested (.ppm)
%%%     - arr{index}.pmra       % Proper Motion RA, mas/y, if requested (.ppm)
%%%     - arr{index}.pmdec      % Proper Motion Dec, mas/y, if requested (.ppm)
%%%     - arr{index}.phot_*     % Light database mean magnitudes, if requested (.mags)
%%%     - arr{index}.source_id  % Source ID from Gaia; 64-bit, if requested (.heavy)
%%%     - arr{index}.*_error    % Standard errors, if requested (.heavy)
%%%     - arr{index}.*_corr     % Correlation coefficients, if requested (.heavy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Sample Usage:
%%%
%%%   octave:1> addpath('jsonlab')     %%% Cf. git clone https://github.com/fangq/jsonlab.git jsonlab/
%%%   octave:2> x.fov={{1,2},.3};      %%% Configure 0.3deg radius circular FOV at RA=1deg, Dec=2deg
%%%   octave:3> x.limit=2;             %%% Configure to return only two lowest (G) magnitude stars
%%%   octave:4> fov_cmd(x)             %%% Get the data ...
%%%
%%%   ans = 
%%%   {
%%%     [1,1] =
%%%
%%%       scalar structure containing the fields:
%%%
%%%         mean_mag =  7.9053
%%%         ra =  1.1257
%%%         dec =  2.2674
%%%         offset =  116655034
%%%
%%%     [1,2] =
%%%
%%%       scalar structure containing the fields:
%%%
%%%         mean_mag =  9.9577
%%%         ra =  1.0022
%%%         dec =  2.2425
%%%         offset =  116655047
%%%
%%%   }
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function return_object = fov_cmd(dotted)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build Python script command

%%% - Required field, .fov
pycmd = 'python fov_cmd.py';
for arg=dotted.fov;
  delim = ' ';
  for val=arg{1,1}
    try
      pycmd=[pycmd delim num2str(val{1,1})];
    catch
      pycmd=[pycmd delim num2str(val(1))];
    end
    delim=',';
  end
end

%%% - Optional fields
try dotted.heavy                                            ;
    pycmd = [ pycmd ' --heavy']                             ; catch end
try dotted.mags                                             ;
    pycmd = [ pycmd ' --mags']                              ; catch end
try dotted.ppm                                              ;
    pycmd = [ pycmd ' --ppm']                               ; catch end
try dotted.j2000                                            ;
    pycmd = [ pycmd ' --j2000']                             ; catch end
try pycmd = [ pycmd ' --limit=' num2str(dotted.limit)]      ; catch end
try pycmd = [ pycmd ' --mag-max=' num2str(dotted.magmax)]   ; catch end
try pycmd = [ pycmd ' --mag-min=' num2str(dotted.magmin)]   ; catch end
try pycmd = [ pycmd ' --mag-type=' dotted.magtype]          ; catch end
try pycmd = [ pycmd ' --gaia-sqlite3=' dotted.gaiasqlite3]  ; catch end
try pycmd = [ pycmd ' --buffer=' num2str(dotted.buffer)]    ; catch end
try pycmd = [ pycmd ' --obsy=' num2str(dotted.obsy)]        ; catch end
try pycmd = [ pycmd ' --obspos=' num2str(dotted.obspos(1))] ;
    pycmd = [ pycmd ','          num2str(dotted.obspos(2))] ;
    pycmd = [ pycmd ','          num2str(dotted.obspos(3))] ; catch end
try pycmd = [ pycmd ' --obsvel=' num2str(dotted.obsvel(1))] ;
    pycmd = [ pycmd ','          num2str(dotted.obsvel(2))] ;
    pycmd = [ pycmd ','          num2str(dotted.obsvel(3))] ; catch end
try pycmd = [ pycmd ' --obsy=' num2str(dotted.obsy)]        ; catch end

[cmd_status,cmd_stdout] = system(pycmd);

fprintf(1,'%s\n',cmd_stdout)

return_object = loadjson(cmd_stdout);

end
