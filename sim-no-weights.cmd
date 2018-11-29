@echo off
REM Generate random model points in interval [0:1,0:1] without initial noise
SET LUA_PATH=lua\?.lua
SET LUA_CPATH=lib/lua/5.1/?.dll
cd %0\..
CALL icpsim.cmd -g "0:1,0:1" -i
