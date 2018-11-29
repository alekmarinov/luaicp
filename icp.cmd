@echo off
REM Starts ICP program
REM SET LUA_PATH=share/lua/?.lua
SET LUA_CPATH=lib/lua/5.1/?.dll
bin\lua5.1.exe lua/icp.lua %*
