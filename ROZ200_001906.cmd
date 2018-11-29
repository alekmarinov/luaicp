@echo off
CD %0\..
CALL icpsim.cmd -i data\ROZ200_001906_normals.txt > data\ROZ200_001906_noised.txt
