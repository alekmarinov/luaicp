-------------------------------------------------------------------------------
-- Implements the Iterative Closest Point algorithm
-- @copyright (C) 2003-2010,  BAS-IIT
-- @author Alexander Marinov (amarinov@iinf.bas.bg)
--         Nadezhda Zlateva  (nzlateva@iinf.bas.bg)
--         Delian Marinov
-------------------------------------------------------------------------------

require "luamatrix"
require "lrun.ai.icp.transform"

local matrix, table, setmetatable, assert, type, icp =
      matrix, table, setmetatable, assert, type, lrun.ai.icp

module "lrun.ai.icp.algorithm"

setmetatable(_M, {__index=matrix})

--- Approaches set A to set B by iterative translations and rotations <br>
-- @param A is the set to be moved
-- @param B is the target set for A movements
-- @param map is a user defined function which returns the subset of B - the closest points to A.
-- @param distance is a user defined function which returns the distance estimation for 2 arbitrary points
-- @param wA is optional point weights if set A
-- @param wB is optional point weights if set B <br>
-- @return iterator function
-- @usage for A, Bs, R, t, wBs in algorithm(A, B, icp.map.one2many, icp.distance.euclid) do ... end <br>
function iterator(A, B, map, distance, wA, wB)
	if not check(A) then A = fromtable(A) end
	if not check(B) then B = fromtable(B) end
	assert(type(map) == "function")
	assert(B:size() >= A:size(), "|B|="..B:size().." must be equal or greater to |A|="..A:size())

	return function ()
		local Bs, wBs = map(A, B, distance, wA, wB)
		local R, t = icp.transform.transform(A, Bs)
		A = icp.transform.move(A, R, t)
		return A, Bs, R, t, wBs
	end
end

return _M
