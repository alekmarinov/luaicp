----------------------------------------------------------------------
--                                                                   --
-- Copyright (C) 2003-2010,  BAS-IIT                                 --
--                                                                   --
-- Project:       Astroinformatics                                   --
-- Filename:      icpsim.lua                                         --
-- Description:   ICP simulation program                             --
-- Authors:       Alexander Marinov (amarinov@iinf.bas.bg)           --
--                Nadezhda Zlateva  (nzlateva@iinf.bas.bg)           --
--                                                                   --
----------------------------------------------------------------------- 

require "luamatrix"
require "lrun.ai.icp.distance"
require "lrun.ai.icp.transform"
require "lrun.ai.icp.map"
require "lrun.ai.icp.algorithm"
local ICP = lrun.ai.icp
local arg = arg

module ("icp.sim", package.seeall)

_NAME = "icpsim"
_VERSION = "1.0"
_DESCRIPTION = "ICP algorithm simulation program"

MAX_ITERS = 40

local appwelcome = _NAME.." ".._VERSION..", ".._DESCRIPTION
local usagetext = "Usage: ".._NAME.." [OPTION]... [INPUTFILE]"
local usagetexthelp = "Try ".._NAME.." --help' for more options."
local errortext = _NAME..": %s"
local helptext = [[
-n   --number N       number of points N to generate if not input file specified
-g   --range  G       scale in range G=X1min:X1max,X2min:X2max,...[,Wmin:Wmax]
-d   --dilute D       diluting deltas D=D1,D2,...[,W], randoms in -Di/2:Di/2
-s   --random-seed S  random seed S, default current time
-t   --translate T    translation vector T=X1,X2,...[,W]
-r   --rotate DEG     rotation DEG degrees
-w   --weighted       use the last column as point weights
-m   --max-iters      maximum iterations limit
-i   --interactive    interactive mode
-pd, --print-distance print the distance at each ICP iteration
-h,  --help           print this help.
]]

POINTSIZE = 5
CROSSSIZE = 0.03
MOVESTEP = 0.01

_opts = nil

--- exit with usage information when the application arguments are wrong 
local function usage(errmsg)
    io.stderr:write("\n"..string.format(errortext, errmsg).."\n")
    io.stderr:write(usagetexthelp)
    os.exit(1)
end

--- exit with error message
local function exiterror(errmsg)
    io.stderr:write(string.format(errortext, errmsg).."\n")
    assert(errmsg)
    os.exit(1)
end

local function parsenums(nums)
	local N = {}
	string.gsub(nums, "([%d%.%-]+)", function (v)
		table.insert(N, tonumber(v))
	end)
	return N
end

--- parses program arguments
local function parseargs(...)
    local opts = {}
    local args = {...}
    local i = 1
    while i <= #args do
	local arg = args[i]
        if arg == "-n" or arg == "--number" then
		opts.npoints = tonumber(args[i+1])
		i = i + 1
        elseif arg == "-g" or arg == "--range" then
		opts.range = {}
		string.gsub(","..args[i+1], "([^,]+)", function(srange)
			local min, max = unpack(parsenums(srange))
			if not (min and max) then
				usage("invalid range format, expected X1min:X1max,X2min:X2max,...[,Wmin:Wmax]")
			end
			table.insert(opts.range, {min=min, max=max})
		end)
		i = i + 1
        elseif arg == "-d" or arg == "--dilute" then
		opts.dilute = parsenums(args[i+1])
		i = i + 1
        elseif arg == "-t" or arg == "--translate" then
		opts.translate = parsenums(args[i+1])
		i = i + 1
        elseif arg == "-r" or arg == "--rotate" then
		opts.rotate = tonumber(args[i+1])
		i = i + 1
        elseif arg == "-s" or arg == "--random-seed" then
		opts.randomseed = tonumber(args[i+1])
		i = i + 1
        elseif arg == "-w" or arg == "--weighted" then
		opts.weighted = true
        elseif arg == "-m" or arg == "--max-iters" then
		opts.maxiters = tonumber(args[i+1])
		i = i + 1
        elseif arg == "-i" or arg == "--interactive" then
		opts.interactive = true
        elseif arg == "-pd" or arg == "--print-distance" then
			opts.print_distance = true
        elseif arg == "-h" or arg == "--help" then
            io.stderr:write(appwelcome.."\n")
            io.stderr:write(usagetext.."\n\n")
            io.stderr:write(helptext)
            os.exit(1)
        else
		if opts.inputfile then
			usage("too many arguments")
		else
			opts.inputfile = arg
		end
        end
	i = i + 1
    end

    if not opts.range then
	if not opts.inputfile then
		opts.range = { {min=0, max=1}, {min=0, max=1} }
	end
    end

    if opts.dilute and #opts.dilute ~= #opts.range then
	usage("diluting deltas must be "..#opts.range)
    end

    if opts.translate and #opts.translate ~= #opts.range then
	usage("translation vector size must be "..#opts.range)
    end

    return opts
end

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

-- return true if the passed values are close to zero
function arezeros(...)
	local eps = matrix.eps
	for _, v in ipairs{...} do
		if math.abs(v) > eps then
			return false
		end
	end
	return true
end

-- return true if rotation matrix is not identity
function isrotate(R)
	return not arezeros(R[1][1]-1, R[2][2]-1, R[1][2], R[2][1])
end

-- return true if translation vector is not identity
function istranslate(t)
	return not arezeros(t[1], t[2])
end

function loadpoints(file, opts)
	local X = {}
	local minmax
	local line = file:read("*l")
	local nline = 1
	while line do
		local point = parsenums(line)
		if opts.range then
			if #opts.range > #point then
				file:close()
				usage("number of point coordinates is less than assigned range dimmension at line #"..nline)
			elseif #opts.range < #point then
				while #point > #opts.range do
					table.remove(point)
				end
			end
		end

		if not minmax then
			minmax = {}
			for i = 1, #point do
				table.insert(minmax, {})
			end
		end

		table.insert(X, point)
		for i = 1, #point do
			local xi = point[i]
			minmax[i].min = minmax[i].min or xi
			minmax[i].max = minmax[i].max or xi
			if xi < minmax[i].min then
				minmax[i].min = xi
			elseif xi > minmax[i].max then
				minmax[i].max = xi
			end
		end
		line = file:read("*l")
		io.stderr:write(string.format("Parsing line %d        \r", nline)) io.stderr:flush()
		nline = nline + 1
		if opts.npoints and nline == opts.npoints then
			break
		end
	end

	opts.range = opts.range or minmax

	-- scale points in opts.range
	for i, point in ipairs(X) do
		for j, coord in ipairs(point) do
			point[j] = opts.range[j].min + (coord - minmax[j].min) * (opts.range[j].max - opts.range[j].min) / (minmax[j].max - minmax[j].min)
		end
	end
	X = matrix.fromtable(X)
	local weights
	if opts.weighted then
		io.stderr:write("Processing data weights...\n") io.stderr:flush()
		local dim = X[1]:size()
		weights = X:col(dim)
		X = #((#X):slice(1, dim-1))
	end
	return X, minmax, weights
end

function savepoints(file, points, weights)
	for i = 1, points:size() do
		for j = 1, points[1]:size() do
			if j > 1 then
				file:write(" ")
			end
			file:write(string.format("%.3f", points[i][j]))
		end
		if weights then
			file:write(string.format(" %.3f", weights[i]))
		end
		file:write("\n")
	end
end

function centroid(X)
	local xc = {}
	for i = 1, X[1]:size() do
		table.insert(xc, X:col(i):sum()/X:size())
	end
	return matrix.fromtable(xc)
end

function distanceset(A, B, distance)
	local res = 0
	for i = 1, A:size() do
		res = res + distance(A[i], B[i])
	end
	return res
end

function rand(min, max)
	if max < min then
		min, max = max, min
	end
	if max - min < matrix.eps then
		return max
	end
	return min + math.random(1000000*(max-min))/1000000
end

-- add diluting noise to X
function noisedilute(X, weights, deltas)
	local M = {}
	local W = {}
	for i = 1, X:size() do
		local p = X[i]
		local tp = {}
		for j = 1, p:size() do
			table.insert(tp, p[j] + rand(-deltas[j]/2, deltas[j]/2))
		end
		table.insert(M, tp)
		if weights then
			table.insert(W, weights[i] + rand(-deltas[X[1]:size() + 1]/2, deltas[X[1]:size() + 1]/2))
		end
	end
	W = weights and matrix.fromtable(W) or nil
	return matrix.fromtable(M), W
end

-- add translational noise to X
function noisetranslate(X, direction)
	if not matrix.check(direction) then direction = matrix.fromtable(direction) end
	local t = {}
	for i = 1, X[1]:size() do
		if direction:size() >= i then
			table.insert(t, direction[i])
		end
	end
	io.stderr:write("translating to "..table.concat(t, ",").."\n")
	local M = {}
	for i = 1, X:size() do
		local p = X:row(i)
		local tp = {}
		for j = 1, p:size() do
			table.insert(tp, p[j] + (t[j] or 0))
		end
		table.insert(M, tp)
	end
	return matrix.fromtable(M)
end

-- add rotational noise to X
function noiserotate(X, degree)
	io.stderr:write("rotating to "..degree.." degree\n")
	alpha = degree * math.pi / 180
	local R = matrix.fromtable{
		{math.cos(alpha), -math.sin(alpha)},
		{math.sin(alpha), math.cos(alpha)},
	}
	local M = {}
	local xc = centroid(X):slice(1, 2)
	for i = 1, X:size() do
		local p = matrix.fromtable{X[i][1], X[i][2]}
		p = R * ( p - xc ) + xc
		local rv = {p[1], p[2]}
		for j = 3, X[1]:size() do
			table.insert(rv, X[i][j])
		end
		table.insert(M, rv)
	end
	return matrix.fromtable(M)
end

-- add random positional noise to set X according to the specified options
function noise(X, weights, opts)
	if opts.randomseed then
		math.randomseed(opts.randomseed)
	else
		math.randomseed(os.time())
	end
	if opts.rotate then
		X = noiserotate(X, opts.rotate)
	end
	if opts.translate then
		X = noisetranslate(X, opts.translate)
	end
	if opts.dilute then
		io.stderr:write("diluting with "..table.concat(opts.dilute, ",").."\n")
		X, weights = noisedilute(X, weights, opts.dilute)
	end
	return X, weights
end

function generatepoints(opts)
	local npoints = opts.npoints or 100
	local X = {}
	local W = {}
	for i = 1, npoints do
		local p = {}
		for j, r in ipairs(opts.range) do
			if not opts.weighted or j < #opts.range then
				table.insert(p, rand(r.min, r.max))
			end
		end
		table.insert(X, p)
		if opts.weighted then
			local weight = rand(opts.range[#opts.range].min, opts.range[#opts.range].max)
			table.insert(W, weight)
		end
	end
	return matrix.fromtable(X), opts.weighted and matrix.fromtable(W)
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

function drawpoints(X, weights, shape, color)
	assert(type(shape) == "string")
	if X then
		if shape == "point" then
			gl.PointSize(POINTSIZE)
			gl.Begin("POINTS")
		elseif shape == "cross" then
			gl.Begin("LINES")
		end
		for i=1, X:size() do
			if weights then
				local a = 1-weights[i]
				color = {a, a, a}
			end
			gl.Color(unpack(color))
			if shape == "point" then
				gl.Vertex(X[i][1], X[i][2])
			elseif shape == "cross" then
				local c = X[i]
				gl.Vertex(c[1]-CROSSSIZE/2, c[2])
				gl.Vertex(c[1]+CROSSSIZE/2, c[2])
				gl.Vertex(c[1], c[2]-CROSSSIZE/2)
				gl.Vertex(c[1], c[2]+CROSSSIZE/2)
			else
				assert(nil, "invalid shape type")
			end
		end
		gl.End()
	end
end

function drawconnectors(X, Y)
	if X and Y then
		gl.Color(0.5, 0.5, 0.5)
		gl.Begin("LINES")
		for i=1, X:size() do
			gl.Vertex(X[i][1], X[i][2])
			gl.Vertex(Y[i][1], Y[i][2])
		end
		gl.End()
	end
end

local points1, weights1, points2, weights2, subpoints, subweights

function _G.DrawFrame()
	gl.MatrixMode("PROJECTION")
	gl.LoadIdentity()
	gl.Ortho(-0.2, 1.2, -0.2, 1.2, -1.0, 1.0)
	gl.MatrixMode("MODELVIEW")
	gl.LoadIdentity()
	gl.ClearColor(1,1,1,1)
	gl.Clear("DEPTH_BUFFER_BIT,COLOR_BUFFER_BIT")
	gl.BlendFunc("SRC_ALPHA", "ONE_MINUS_SRC_ALPHA")
	gl.Enable("BLEND")
	gl.Enable("LINE_SMOOTH") 
	
	drawpoints(points1, weights1, "point", {0, 0, 0})
	drawpoints(points2, weights2, "cross", {0, 0, 0})

	if subpoints then
		drawconnectors(points2, subpoints)
	end

	glut.SwapBuffers()
	gl.Flush()
end

function _G.Reshape(width, height)

	gl.Viewport(0, 0, width, height)

	gl.MatrixMode('PROJECTION')
	gl.LoadIdentity()

	gl.MatrixMode('MODELVIEW')
	gl.LoadIdentity()

	DrawFrame()
end

function _G.Keyboard(key)
	if key == 27 then
		os.exit()
	elseif key == string.byte("d") then
		local size = points2:size()
		local idx = math.random(size-1)
		points2 = points2:slice(1, idx-1)..points2:slice(idx+1, size)
		if weights2 then
			local part1 = weights2:slice(1, idx-1)
			local part2 = weights2:slice(idx+1, size)
			local W = {}
			for i = 1, part1:size() do table.insert(W, part1[i]) end
			for i = 1, part2:size() do table.insert(W, part2[i]) end

			weights2 = matrix.fromtable(W)
		end
		subpoints, subweights = ICP.map.one2many(points2, points1, ICP.distance.weighted, weights2, weights1)
		io.stderr:write("size reduced to "..points2:size().."\n")
		io.stderr:flush()
	elseif key == string.byte("q") then
		points2 = noiserotate(points2, 1, true)
	elseif key == string.byte("w") then
		points2 = noiserotate(points2, -1, true)
	elseif key == string.byte("j") then
		points2 = noisetranslate(points2, matrix.fromtable{-MOVESTEP, 0}, true)
	elseif key == string.byte("i") then
		points2 = noisetranslate(points2, matrix.fromtable{0, MOVESTEP}, true)
	elseif key == string.byte("l") then
		points2 = noisetranslate(points2, matrix.fromtable{MOVESTEP, 0}, true)
	elseif key == string.byte("k") then
		points2 = noisetranslate(points2, matrix.fromtable{0, -MOVESTEP}, true)
	elseif key == string.byte("m") then
		subpoints, subweights = ICP.map.one2many(points2, points1, ICP.distance.weighted, weights2, weights1)
	elseif key == string.byte("o") then
		io.stderr:write("One 2 one mapping started\n") io.stderr:flush()
		subpoints, subweights = ICP.map.one2one(points2, points1, ICP.distance.weighted, weights2, weights1)
		io.stderr:write("One 2 one mapping finished\n") io.stderr:flush()
	elseif key == string.byte("t") then
		if subpoints then
			local R, t = ICP.transform(points2, subpoints, weights2, subweights)
			points2 = ICP.move(points2, R, t)
		end
	elseif key == string.byte("r") then
		local deltas = {0.01, 0.01}
		if weights2 then
			table.insert(deltas, 0.01)
		end
		points2, weights2 = noisedilute(points2, weights2, deltas)
	elseif key == 32 then
		local i = 1
		local delay1 = false
		-- Note: local minimum is detected if no translation nor rotation occurred at 2 iterations
		--       because of possible new mapping after the 1st one
		io.stderr:write(string.format("Performing ICP on %d data over %d model points...\n", points2:size(), points1:size())) io.stderr:flush()
		if _opts.print_distance then print(distanceset(points2, points1, ICP.distance.weighted)) end
		local maxiters = _opts.maxiters or MAX_ITERS
		for data, model, R, t, modelweights in ICP.algorithm.iterator(points2, points1, ICP.map.one2many, ICP.distance.weighted, weights2, weights1) do
			points2 = data
			if _opts.print_distance then print(distanceset(points2, points1, ICP.distance.weighted)) end
			subpoints = model
			subweights = modelweights
			DrawFrame()
			if i == maxiters then
				print("Stop on exceeding maximum ("..maxiters..") iterations")
				break
			elseif delay1 then
				print("Stops on local minimum in "..i.." iterations")
				break
			elseif not (istranslate(t) or isrotate(R)) then
				delay1 = true
			end
			i = i + 1
		end
	elseif key == string.byte("s") then
		savepoints(io.stdout, points2, weights2)
		io.stdout:flush()
	elseif key == 13 then
		points2 = points1:copy()
		if weights then
			weights2 = weights1:copy()
		end
	elseif key == string.byte("h") then
		io.stderr:write([[
d - remove random point
q - rotate right
w - rotate left
i - move up
j - move right
k - move left
l - move down
m - map one to many
o - map one to one
t - transform
r - dilute points randomly with 0.01 offsets (including weights)
s - write set to stdout
space - ICP 50 times
enter - reset state
esc - quit
]])
		io.stderr:flush()
	end
	DrawFrame()
end

function InitGL()
	require "gl"

	glut.Init()
	glut.InitDisplayMode("RGB,SINGLE")
	glut.InitWindowSize(600, 600)
	glut.CreateWindow(_NAME.." ".._VERSION)
	glut.DisplayFunc('DrawFrame')
	glut.KeyboardFunc('Keyboard')
	glut.ReshapeFunc('Reshape')
	--glut.FullScreen(true)
end

function main(...)
	--- parse program arguments
	local opts = parseargs(...)
	_opts = opts

	local points, minmax, weights
	if opts.inputfile then
		local file = assert(io.open(opts.inputfile, "r"))
		points, minmax, weights = loadpoints(file, opts)
		io.stderr:write(opts.inputfile.." loaded\n") io.stderr:flush()
		file:close()
		points = normalize(points, minmax)
	else
		points, weights = generatepoints(opts)
	end

	if opts.interactive then
		points1 = points
		weights1 = weights
		points2 = points:copy()
		if weights then
			weights2 = weights:copy()
		end
		points2, weights2 = noise(points2, weights2, opts)
		InitGL()
		glut.PostRedisplay()
		glut.MainLoop()
	else
		savepoints(io.stdout, points, weights)

		if minmax then
			io.stderr:write("input range = ")
			for i, mm in ipairs(minmax) do
				if i > 1 then
					io.stderr:write(",")
				end
				io.stderr:write(mm.min..":"..mm.max)
			end
			io.stderr:write("\n")
			io.stderr:flush()
		end
	end
end

main(unpack(arg))

return _M
