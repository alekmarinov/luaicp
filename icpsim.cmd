@echo off
REM Starts ICP simulation program
SET LUA_PATH=lua\?.lua
SET LUA_CPATH=lib/lua/5.1/?.dll
cd %0\..
bin\lua5.1.exe lua/icpsim.lua %*
