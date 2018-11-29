@echo off
CD %0\..
CALL icpsim.cmd -pd -s 1234 -t 0.2,0.2,0 -w -r 15 -g 0:1,0:1,0:1 -i data\ROZ200_001906_mag.txt
