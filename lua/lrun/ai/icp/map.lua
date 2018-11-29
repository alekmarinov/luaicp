-------------------------------------------------------------------------------
-- Implements several mapping functions between 2 sets
-- @copyright (C) 2003-2010,  BAS-IIT
-- @author Alexander Marinov (amarinov@iinf.bas.bg)
--         Nadezhda Zlateva  (nzlateva@iinf.bas.bg)
--         Delian Marinov
-------------------------------------------------------------------------------

require "luamatrix"

local matrix, table, setmetatable, assert, type, print, icp =
      matrix, table, setmetatable, assert, type, print, lrun.ai.icp

module "lrun.ai.icp.map"

setmetatable(_M, {__index=matrix})

--- Computes subset of B - the closest points to A
-- such as each element of A has exactly one unique point in B
-- @param A is the set to be mapped to
-- @param B is the set to be mapped from
-- @param distance is a distance function
-- @param wA optional scalar weights of the points in A
-- @param wB optional scalar weights of the points in B
-- @return subset of B
-- @return subset of wB or nil
-- @see lrun.ai.icp.distance
function one2one(A, B, distance, wA, wB)
	if not check(A) then A = fromtable(A) end
	if not check(B) then B = fromtable(B) end
	assert(B:size() >= A:size(), "|B|="..B:size().." must be equal or greater to |A|="..A:size())
	distance = distance or icp.distance.euclid

	local function closestpoints(X1, X2, Ha, Hb)
		local c1, c2
		local d = 9999999
		for i = 1, X1:size() do
			if not Ha[i] then
				for j = 1, X2:size() do
					if not Hb[j] then
						local dij = distance(X1[i], X2[j], wA and wA[i], wB and wB[j])
						if dij < d then
							d = dij
							c1 = i
							c2 = j
						end
					end
				end
			end
		end
		return c1, c2
	end

	local Ha, Hb = {}, {}
	local C = {}
	local W = {}
	for k = 1, A:size() do
		table.insert(C, {})
		if wB then
			table.insert(W, wB[k])
		end
	end
	for k = 1, A:size() do
		local ia, ib = closestpoints(A, B, Ha, Hb)
		Ha[ia] = true
		Hb[ib] = true

		C[ia] = {}
		for j = 1, B[1]:size() do
			table.insert(C[ia], B[ib][j])
		end
		if wB then
			W[ia] = wB[ib]
		end
	end
	return fromtable(C), wB and fromtable(W)
end

-- returns subset of B - the closest points to A
-- one to many mapping variant

--- Computes subset of B - the closest points to A
-- such as each element of A can have one or more corresponding points in B
-- @param A is the set to be mapped to
-- @param B is the set to be mapped from
-- @param distance is a distance function
-- @param wA optional scalar weights of the points in A
-- @param wB optional scalar weights of the points in B
-- @return subset of B
-- @return subset of wB or nil
-- @see lrun.ai.icp.distance
function one2many(A, B, distance, wA, wB)
	if not check(A) then A = fromtable(A) end
	if not check(B) then B = fromtable(B) end
	assert(B:size() >= A:size(), "|B|="..B:size().." must be equal or greater to |A|="..A:size())
	distance = distance or icp.distance.euclid
	local C = {}
	local W = {}
	for i = 1, A:size() do
		local d = distance(A[i], B[i], wA and wA[i], wB and wB[i])
		local closest = i
		for j = 1, B:size() do
			local dij = distance(A[i], B[j], wA and wA[i], wB and wB[j])
			if dij < d then
				d = dij
				closest = j
			end
		end
		local v = {}
		for j = 1, B[1]:size() do
			table.insert(v, B[closest][j])
		end
		table.insert(C, v)
		if wB then
			table.insert(W, wB[closest])
		end
	end
	return fromtable(C), wB and fromtable(W)
end

return _M
