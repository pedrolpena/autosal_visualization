addpath(genpath([pwd(),filesep,'m']));

if is_octave ==1
    warning ('off', 'Octave:divide-by-zero');
    graphics_toolkit('gnuplot');
    set (0, 'defaultaxesfontname', 'Helvetica');
    pkg load statistics
    more off
end