----------------------------------------------------------------------
--                                                                   --
-- Copyright (C) 2003-2010,  BAS-IIT                                 --
--                                                                   --
-- Project:       Astroinformatics                                   --
-- Filename:      icp.lua                                            --
-- Description:   Align two point sets with ICP algorithm            --
-- Authors:       Alexander Marinov (amarinov@iinf.bas.bg)           --
--                Nadezhda Zlateva  (nzlateva@iinf.bas.bg)           --
--                                                                   --
----------------------------------------------------------------------- 

require "gl"
require "luamatrix"
require "lrun.ai.icp.distance"
require "lrun.ai.icp.transform"
require "lrun.ai.icp.map"
require "lrun.ai.icp.algorithm"
local ICP = lrun.ai.icp
local arg = arg

module ("icp", package.seeall)

_NAME = "icp"
_VERSION = "1.0"
_DESCRIPTION = "align two point sets with ICP algorithm"

local appwelcome = _NAME.." ".._VERSION..", ".._DESCRIPTION
local usagetext = "Usage: ".._NAME.." [OPTION]... [MODELFILE] [DATAFILE]"
local usagetexthelp = "Try ".._NAME.." --help' for more options."
local errortext = _NAME..": %s"
local helptext = [[
-h,  --help                      print this help.
]]

State = {
	pointsize = 5,
	crosssize = 0.03,
	scale = 1
}

-- print matrix or vector
local function mprint(m, l)
	if l then print (l) end
	  assert(matrix.check(m), 'matrix expected, got ' .. type(m))
	  assert(matrix.dim(m) <= 2)
	  -- maximum entry length
	  local ml = matrix.reduce(m, function(x, y)
		local s = tostring(y):len();
		if x < s then x = s end
		return x
	  end, 0)
	  if matrix.dim(m) == 1 then -- vector?
		matrix.map(m, function(e)
			local es = tostring(e)
			print(string.rep(' ', 2 + ml - es:len()) .. e)
		end)
	  else  -- array
		for i = 1, m:size(1) do -- print each row
			print(matrix.reduce(m[i], function(x, y)
				local es = tostring(y)
				return x .. string.rep(' ', 2 + ml - es:len()) .. es
			end, ''))
		end
	  end
end

--- exit with usage information when the application arguments are wrong 
local function usage(errmsg)
    io.stderr:write(string.format(usagetext, errmsg).."\n\n")
    io.stderr:write(usagetexthelp)
    os.exit(1)
end

--- exit with error message
local function exiterror(errmsg)
    io.stderr:write(string.format(errortext, errmsg).."\n")
    assert(errmsg)
    os.exit(1)
end

--- parses program arguments
local function parseargs(...)
    local opts = {}
    local args = {...}
    for i, arg in ipairs(args) do
        if arg == "-h" or arg == "--help" then
            io.stderr:write(appwelcome.."\n")
            io.stderr:write(usagetext.."\n\n")
            io.stderr:write(helptext)
            os.exit(1)
        else
		if not opts.filemodel then
			opts.filemodel = arg
		elseif not opts.filedata then
			opts.filedata = arg
		else
			usage("too many arguments")
		end
        end
    end
    if not opts.filemodel then
        usage("missing MODELFILE")
    elseif not opts.filedata then
        usage("missing DATAFILE")
    end

    return opts
end

function normalize(X)
	local M = {}
	local mins, maxs = {}, {}
	for i = 1, X[1]:size() do
		local min, max = X:col(i):min(), X:col(i):max()
		table.insert(mins, min)
		table.insert(maxs, max)
	end
	for i = 1, X:size() do
		local p = {}
		for j = 1, X[1]:size() do
			table.insert(p, (X[i][j] - mins[j])/(maxs[j] - mins[j]))
		end
		table.insert(M, p)
	end
	return matrix.fromtable(M)
end

function loadpoints(filename)
	local file, err = io.open(filename)
	if not file then
		return nil, err
	end
	local function parsenums(nums)
		local N = {}
		string.gsub(nums, "([%d%.%-]+)", function (v)
			table.insert(N, tonumber(v))
		end)
		return N
	end

	local X = {}
	local line = file:read("*l")
	while line do
		local point = parsenums(line)
		table.insert(X, point)
		line = file:read("*l")
	end
	file:close()

	return normalize(matrix.fromtable(X))
end

function centroid(X)
	return matrix.fromtable{X:col(1):sum(), X:col(2):sum()}/X:size()
end

function drawcentroid(c, color, shape)
	if c then
		shape = shape or "point"
		gl.Color(color)
		if shape == "point" then
			gl.PointSize(State.pointsize*2)
			gl.Begin("POINTS")
			gl.Vertex(c[1], c[2])
			gl.End()
			gl.PointSize(State.pointsize)
			gl.Color(1, 1, 1)
			gl.Begin("POINTS")
			gl.Vertex(c[1], c[2])
			gl.End()
		elseif shape == "cross" then
			gl.Begin("LINES")
			gl.Vertex(c[1]-State.crosssize, c[2])
			gl.Vertex(c[1]+State.crosssize, c[2])
			gl.Vertex(c[1], c[2]-State.crosssize)
			gl.Vertex(c[1], c[2]+State.crosssize)
			gl.End()
		end
	end
end

function drawpoints(X, shape, color)
	assert(type(shape) == "string")
	if X then
		if shape == "point" then
			gl.PointSize(State.pointsize)
			gl.Begin("POINTS")
		elseif shape == "cross" then
			gl.Begin("LINES")
		end
		for i=1, X:size() do
			if X[1]:size() == 3 then
				local a = X[i][3]
				color = {a, a, a}
			end
			gl.Color(unpack(color))
			if shape == "point" then
				gl.Vertex(X[i][1], X[i][2])
			elseif shape == "cross" then
				local c = X[i]
				gl.Vertex(c[1]-State.crosssize/2, c[2])
				gl.Vertex(c[1]+State.crosssize/2, c[2])
				gl.Vertex(c[1], c[2]-State.crosssize/2)
				gl.Vertex(c[1], c[2]+State.crosssize/2)
			else
				assert(nil, "invalid shape type")
			end
		end
		gl.End()
	end
end

function drawconnectors(X, Y, color)
	if X and Y then
		gl.Color(color)
		gl.Begin("LINES")
		for i=1, X:size() do
			gl.Vertex(X[i][1], X[i][2])
			gl.Vertex(Y[i][1], Y[i][2])
		end
		gl.End()
	end
end

function _G.Reshape(width, height)

	gl.Viewport(0, 0, width, height)

	gl.MatrixMode('PROJECTION')
	gl.LoadIdentity()

	gl.MatrixMode('MODELVIEW')
	gl.LoadIdentity()

	DrawFrame()
end

function _G.DrawFrame()
	gl.MatrixMode("PROJECTION")
	gl.LoadIdentity()
	--local cx, cy = (State.range.maxx - State.range.minx) / 20, (State.range.maxy - State.range.miny) / 20
	--gl.Ortho(State.range.minx-cx, State.scale*State.range.maxx+cx, State.scale*State.range.maxy+cy, State.scale*State.range.miny-cy, -1.0, 1.0)
	gl.Ortho(-0.2, 1.2, -0.2, 1.2, -1.0, 1.0)
	gl.MatrixMode("MODELVIEW")
	gl.LoadIdentity()
	gl.ClearColor(1,1,1,1)
	gl.Clear("DEPTH_BUFFER_BIT,COLOR_BUFFER_BIT")
	gl.BlendFunc("SRC_ALPHA", "ONE_MINUS_SRC_ALPHA")
	gl.Enable("BLEND")
	gl.Enable("LINE_SMOOTH") 
	
	-- draw model
	drawpoints(State.model, "point", {0.7, 0.7, 0.7})
	--drawpoints(State.submodel, "point", {0.3, 0.3, 0.3})
	--drawcentroid(State.submodelcenter, {0, 0, 0, 0.7})
	
	-- draw data
	drawpoints(State.data, "cross", {0, 0, 0})
	--drawcentroid(State.datacenter, "cross", {0, 0, 0})

	-- connectors of data with submodel
	drawconnectors(State.data, State.submodel, {0, 0, 0, 0.7})

	glut.SwapBuffers()
	gl.Flush()
end

_G.Idle = DrawFrame

function Run()
	local i = 1
	for Data, Model, Ri, ti in ICP.algorithm(State.data, State.model, ICP.map.one2many) do
		print(i..":"..ICP.distance.euclid(Data, Model))
		State.data = Data
		State.submodel = Model
		State.datacenter = centroid(Data)
		State.submodelcenter = centroid(Model)
		i = i + 1
		DrawFrame()
		if i == 200 then
			break
		end
	end
end

function _G.Keyboard(key)
	if key == 27 then os.exit() end
	if key == 49 then -- 1
		State.submodel = ICP.map.one2many(State.data, State.model)
		State.submodelcenter = centroid(State.submodel)
		State.datacenter = centroid(State.data)
		print("findmatch d = "..ICP.distance(State.submodel, State.data))
		glut.PostRedisplay()
	end
	if key == 50 then -- 2
		local R, t = ICP.transform(State.data, State.submodel)
		mprint(R, "R")
		mprint(t, "t")
		State.data = ICP.move(State.data, R, t)
		State.datacenter = centroid(State.data)
		print("findtransform d = "..ICP.distance(State.submodel, State.data))
		glut.PostRedisplay()
	end
	if key == 32 then
		Run()
	end
end

function InitGL()
	glut.Init()
	glut.InitDisplayMode("RGB,SINGLE")
	glut.InitWindowSize(600, 600)
	glut.CreateWindow(_NAME.." ".._VERSION)
	glut.DisplayFunc('DrawFrame')
	glut.KeyboardFunc('Keyboard')
	glut.IdleFunc('Idle')
	glut.ReshapeFunc('Reshape')
	--glut.FullScreen(true)
end

function main(...)
	--- parse program arguments
	local opts = parseargs(...)

	State.model = assert(loadpoints(opts.filemodel))
	State.data = assert(loadpoints(opts.filedata))	

	mprint(State.model, "Model")
	mprint(State.data, "Data")

	print(string.format("loaded %d model points and %d data points", State.model:size(), State.data:size()))
	--print("Press '1' to match, '2' to transform and Space to start algorithm")

	InitGL()
	-- Run()
	glut.MainLoop()
end

main(unpack(arg))

return _M
