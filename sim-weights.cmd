@echo off
REM Generate random weighted model points in interval [0:1,0:1,0:1] and dilute data points with offsets [0:1,0:1,0:0]
SET LUA_PATH=lua\?.lua
SET LUA_CPATH=lib/lua/5.1/?.dll
CD %0\..
CALL icpsim.cmd -d "0.01,0.01,0" -t "0.5,0.5,0" -r 180 -w -g "0:1,0:1,0:1" -i
