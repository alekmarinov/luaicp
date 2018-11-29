-------------------------------------------------------------------------------
-- Implements linear transformation functions
-- @copyright (C) 2003-2010,  BAS-IIT
-- @author Alexander Marinov (amarinov@iinf.bas.bg)
--         Nadezhda Zlateva  (nzlateva@iinf.bas.bg)
--         Delian Marinov
-------------------------------------------------------------------------------

require "luamatrix"

local matrix, table, setmetatable, assert =
      matrix, table, setmetatable, assert

module "lrun.ai.icp.transform"

setmetatable(_M, {__index=matrix})

--- Computes rotation R and translation t aligning A to B
-- such that min|(R*A-t) - B|
-- @param A is the source set to align from
-- @param B is the target set to align to <br>
-- @return rotation R
-- @return translation t
-- @return details table with internaly computed results, e.g. the centroids of A and B 
function transform(A, B)
	if not check(A) then A = fromtable(A) end
	if not check(B) then B = fromtable(B) end
	assert(A[1]:size() == B[1]:size(), "dim(A)="..A[1]:size().." must be equal to dim(B)="..B[1]:size())
	assert(B:size() >= A:size(), "|B|="..B:size().." must be equal or greater to |A|="..A:size())

	local n = A:size()

	-- ac, bc - centroids
	local ac = centroid(A)
	local bc = centroid(B)

	-- a1 = a - bc, b1 = b - ac, a(-B, b(-A
	-- H = Sum(a1[i]*#b1[i])
	local H = zeros(B[1]:size())
	for i = 1, n do
		local a0 = A[i] - ac
		local b0 = B[i] - bc

		-- represent vectors as horizontal matrices
		a1 = fromtable{{a0[1]}}
		b1 = fromtable{{b0[1]}}
		for i = 2, a0:size() do
			a1 = a1..a0:slice(i)
			b1 = b1..b0:slice(i)
		end

		H = H + a1 * #b1
	end

	-- compute the singular value decomposition of H
	local U, S, V = svd(H)

    -- compute R such that maximizes Trace(RH)
	local R = V*#U

	-- compute T minimizing square root difference of A, B
	local t = R * ac - bc

	local details = {
		ac = ac,
		bc = bc,
		H = H,
		U = U,
		S = S,
		V = V
	}

	return R, t, details
end

--- Moves set X by rotation R and translation t
-- @param X is the set to be moved
-- @param R is rotation matrix
-- @param t is translation vector
-- @return new set representing X rotated by R and translated by t 
function move(X, R, t)
	if not check(X) then X = fromtable(X) end
	local M = {}
	for i = 1, X:size() do
		local p = R * X[i] - t
		local v = {}
		for i = 1, p:size() do
			table.insert(v, p[i])
		end
		table.insert(M, v)
	end
	return fromtable(M)
end

--- Computes the centroid of set X
-- @param X is the input set
-- @return vector representing the centroid of X
function centroid(X)
	if not check(X) then X = fromtable(X) end
	local xc = {}
	for i = 1, X[1]:size() do
		table.insert(xc, X:col(i):sum()/X:size())
	end
	return matrix.fromtable(xc)
end

--- Computes the mass center of set X
-- @param X is the input set
-- @param wx is a vector with weights of the points in X
-- @return vector representing the mass center of X
function masscenter(X, wx)
	if not check(X) then X = fromtable(X) end
	assert(X:size() == wx:size())
	local xc = {}
	for i = 1, X[1]:size() do
		table.insert(xc, wx * X:col(i))
	end
	return matrix.fromtable(xc) / wx:sum()
end

return _M
