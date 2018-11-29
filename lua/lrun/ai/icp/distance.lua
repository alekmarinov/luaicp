-------------------------------------------------------------------------------
-- Implements several distance functions
-- @copyright (C) 2003-2010,  BAS-IIT
-- @author Alexander Marinov (amarinov@iinf.bas.bg)
--         Nadezhda Zlateva  (nzlateva@iinf.bas.bg)
--         Delian Marinov
-------------------------------------------------------------------------------

require "luamatrix"

local matrix, math, setmetatable, assert, print =
      matrix, math, setmetatable, assert, print

module "lrun.ai.icp.distance"

setmetatable(_M, {__index=matrix})

--- Computes the square distance between vector a and vector b
-- @param a is the 1st vector
-- @param b is the 2nd vector such as |a| = |b|
-- @return Sum (a[i]-b[i])^2, i = 1..|a|
function square(a, b)
	if not check(a) then a = fromtable(a) end
	if not check(b) then b = fromtable(b) end
	assert(a:size() == b:size(), "expected |a|="..a:size().." == |b|="..b:size())

	local c = a - b

	local d = 0
	for i = 1, c:size() do
		d = d + c[i]*c[i]
	end
	return d
end

--- Computes the euclid distance between vector a and vector b
-- @param a is the 1st vector
-- @param b is the 2nd vector such as |a| = |b|
-- @return ||a-b||
function euclid(a, b)
	return norm(a - b)
end

--- Computes a weighted distance between vector a and vector b<br>
-- if input weights are omited the result distance is the same as euclid 
-- @param a is the 1st vector
-- @param b is the 2nd vector such as |a| = |b|
-- @param wa is optional scalar weight for vector a
-- @param wb is optional scalar weight for vector b
-- @return (Max(wa, wb)/Min(wa, wb))^2 * ||a-b||
function weighted(a, b, wa, wb)
	if not check(a) then a = fromtable(a) end
	if not check(b) then b = fromtable(b) end
	assert(a:dim() == 1, "expected dim(a)=="..1)
	assert(b:dim() == 1, "expected dim(b)=="..1)
	assert(a:size() == b:size(), "expected |a|="..a:size().." == |b|="..b:size())

	local wab = 1
	if wa and wb then
		wab = math.max(wa, wb)/math.min(wa, wb)
	end
	return wab^2*norm(a-b)
end
