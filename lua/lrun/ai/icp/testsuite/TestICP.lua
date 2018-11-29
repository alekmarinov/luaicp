-----------------------------------------------------------------------
--                                                                   --
-- Copyright (C) 2003-2010,  BAS-IIT                                 --
--                                                                   --
-- Project:       Astroinformatics                                   --
-- Filename:      TestICP.lua                                        --
-- Description:   Tests Iterative Closest Point algorithm            --
-- Authors:       Alexander Marinov (amarinov@iinf.bas.bg)           --
--                Nadezhda Zlateva  (nzlateva@iinf.bas.bg)           --
--                                                                   --
----------------------------------------------------------------------- 

require "luamatrix"

require "luaunit"
require "lrun.ai.icp.distance"
require "lrun.ai.icp.transform"
require "lrun.ai.icp.map"
require "lrun.ai.icp.algorithm"
local icp = lrun.ai.icp

module ("TestICP", package.seeall)

local A = matrix.fromtable{
	{-1, 0},
	{0, 0},
	{0, 2},
}

local wA = matrix.fromtable{
	0.1,
	0.5,
	0.9
}

local B = matrix.fromtable{
	{1, 1},
	{1, 0},
	{3, 0},
}


local wB = matrix.fromtable{
	0.1,
	0.5,
	0.9
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

function _testTransform()
	local R, t, details = icp.transform.transform(A, B)
	assertEquals(-1/3, details.ac[1])
	assertEquals(2/3, details.ac[2])
	assertEquals(5/3, details.bc[1])
	assertEquals(1/3, details.bc[2])
	local U, S, V, H = details.U, details.S, details.V, details.H
	local EPS = matrix.zeros(H:size()):fill(matrix.eps)
	local E = matrix.eye(H:size())
	assertTrue((U*S*V - H) <= EPS)
	assertTrue((U*#U - E) <= EPS)
	assertTrue((V*#V - E) <= EPS)
	assertTrue((R*#R - E) <= EPS)
	assertTrue(matrix.isdiagonal(S))
	local A1 = icp.transform.move(A, R, t)
	for i = 1, A1:size() do
		for j = 1, A1[1]:size() do
			assertTrue(math.abs(B[i][j] - A1[i][j]) <= 2*matrix.eps)
		end
	end
	assertTrue(icp.distance.square(A1, B) <= 2*matrix.eps)
end

function _testMap()
	local EPS = matrix.fromtable{matrix.eps, matrix.eps}

	local B1 = icp.map.one2many(A, B)
	assertTrue((B1[1] - B[2]) <= EPS)
	assertTrue((B1[2] - B[2]) <= EPS)
	assertTrue((B1[3] - B[1]) <= EPS)

	local B1 = icp.map.one2one(A, B)
	assertTrue((B1[2] - B[2]) <= EPS)
	assertTrue((B1[3] - B[1]) <= EPS)
	assertTrue((B1[1] - B[3]) <= EPS)
end

function testDistance()
	local A1 = icp.transform.move(A, matrix.eye(2), matrix.fromtable{2, 0})
	assertEquals(12, icp.distance.square(A, A1))

	for i = 1, A:size() do
		for j = 1, B:size() do
			local a, b, wa, wb = A[i], B[j], wA[i], wB[j]
			assertTrue(icp.distance.weighted(a, b, wa, wb) >= icp.distance.euclid(a, b, wa, wb))
		end
	end
end

function _testAlgorithm()
	math.randomseed(1234)
	local A = {}
	local B = {}
	local Range = 50
	for i = 1, 100 do
		local x, y = math.random(Range), math.random(Range)
		table.insert(A, {x, y})
		table.insert(B, {x, y})
	end
	local alpha = math.pi*15/180
	B = icp.transform.move(B,
		matrix.fromtable{ {math.cos(alpha), -math.sin(alpha)}, {math.sin(alpha), math.cos(alpha)} },
		matrix.fromtable{3, 3}
	)

	local dists = {0, 0}
	for i, mapper in ipairs{icp.map.one2one, icp.map.one2many} do
		local iters = 100

		for A1, B1, R1, t1 in icp.algorithm.icp(A, B, mapper) do
			local dist = icp.distance.square(A1, B1)
			io.write(string.format("distance = %.5f           \r", dist))
			if dist == dists[i] or dist < matrix.eps then
				dists[i] = dist
				break
			end
			dists[i] = dist
			iters = iters - 1
			if iters == 0 then break end
		end
	end
	assertTrue(dists[1] <= matrix.eps)
	assertTrue(dists[2] <= matrix.eps)
end

LuaUnit:run("TestICP")
