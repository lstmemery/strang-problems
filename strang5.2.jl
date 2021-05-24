### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 8993d71d-9c83-4b2b-87be-fe5afc2fd359
using Symbolics

# ╔═╡ 45411a2e-449a-4475-838b-187adabfc6b1
using LinearAlgebra

# ╔═╡ 60bb7f84-877a-4e32-b50e-8c63a5c2b24b
using MatrixDepot

# ╔═╡ 16717b6f-4e60-41d6-b796-f2b1c09decc8
using Permutations

# ╔═╡ aa4a53fa-7168-408b-92ed-6c6bb4a4d128
using LaTeXStrings

# ╔═╡ 48f45c84-2fae-4c89-9d1d-b5fdff136ffd
using Markdown

# ╔═╡ 718e651f-a9e6-4e6d-ab50-9a429179bdf8
using Latexify

# ╔═╡ 3b1025b8-ab31-41d2-9ff5-1d31e7e1ecf2
using Memoization

# ╔═╡ 111fa6a0-3d84-420e-b546-59075a03a5c2
begin
	@variables x a b c d e f g h i j k l m n o p q r s t u v w y z
end;

# ╔═╡ e4c74c4c-97c4-4865-bc54-cddf00fa2e58
begin
	@variables a12 a13 a14 a21 a23 a24 a31 a32 a34 a41 a42 a43
end;

# ╔═╡ c3e4419b-7421-423f-a804-0b98f7638598
begin
	@variables a1 a2 b1 b2
end;

# ╔═╡ 742113f5-7cca-4139-b288-a0dc41264423
show_num(i) = L"%$(latexify(i))"

# ╔═╡ 5b4a23b6-17ee-4ce9-814f-cb8dd2489a58
rat(a) = rationalize(a, tol=1e-4)

# ╔═╡ 49f2380f-af3b-4b20-afde-1b793499cdec
md"""
``\newcommand\m[1]{\begin{bmatrix}#1\end{bmatrix}}``
"""

# ╔═╡ e9ed900c-a286-4e29-adaf-5fff401068f9
md"""
``\newcommand\v[1]{\begin{vmatrix}#1\end{vmatrix}}``
"""

# ╔═╡ f880b525-d950-472a-a966-9313c3ad9e14
function big_determinant_terms(A::Array)
	matrix_size = size(A)[1]
	little_det(perm::Permutation, A::Array) = sign(perm) * prod(diag(Matrix(perm)*A))
	return [little_det(Permutation(matrix_size, i), A) for i=1:factorial(matrix_size)]
end

# ╔═╡ 95074837-fbab-43ba-812c-680895172946
big_determinant(A::Array) = sum(big_determinant_terms(A))

# ╔═╡ e210ddda-0ebc-42fe-88b5-6545838d54e1
cofactor_matrix(A, i, j) = A[1:end .!= i, 1:end .!= j]

# ╔═╡ fa7f4e03-f93e-401a-9280-0f94957d22ac
function cofactor_determinant(A)
	if size(A)[1] == 1
		return A[1, 1]
	else
		return sum([A[1,j]*(((-1)^(1 + j)) * cofactor_determinant(cofactor_matrix(A,1,j))) for j=1:size(A)[2]])
	end
end

# ╔═╡ 400ee06e-1172-4805-abf9-981a2a66dadf
cofactors(A) = [(((-1)^(1 + j)) * cofactor_determinant(cofactor_matrix(A,1,j))) for j=1:size(A)[2]]

# ╔═╡ b93f597c-fae2-4ed2-ae24-b2f17a94d8c8
all_cofactors(A) = [(((-1)^(i + j)) * cofactor_determinant(cofactor_matrix(A,i,j))) for i=1:size(A)[1], j=1:size(A)[2]]

# ╔═╡ 7ae07d75-d7e4-467c-8a0a-da084de0d5ac
A = [
	1 2;
	3 4
]

# ╔═╡ aa4721fa-046c-4d28-a8f6-245d4e73a7d0
all_cofactors(A)

# ╔═╡ 862e32ee-e26a-4cd5-ae0f-83a24b898102
big_determinant(A)

# ╔═╡ 1b83dbef-8bc2-4fe0-a357-7f95b9a2d455
cofactor_determinant(A)

# ╔═╡ 9d4a50b6-8fe5-4611-aafe-af548bc5970d
cofactors(A)

# ╔═╡ 95050319-88ab-4c0c-9406-d1dcc36c2ccd
rank(A)

# ╔═╡ 2b5d4b30-aeab-11eb-0b6d-5521705ce8f5
md"""
## Q1
"""

# ╔═╡ b18b009f-f2ea-4093-a91c-1bf56c661a86
q1A = [
	1 2 3;
	3 1 2;
	3 2 1
]

# ╔═╡ 68fe8dde-2fb3-46b6-9d4f-cc9f85c78ee2
big_determinant(q1A)

# ╔═╡ 642bef99-6b5f-41c3-84c0-98a47ab53f21
cofactor_determinant(q1A)

# ╔═╡ 6cce2d8f-2d6d-46e4-917b-a06078acc12a
det(q1A)

# ╔═╡ 33714c08-90c6-4044-ad80-c5649024f47f
q1B = [
	1 2 3;
	4 4 4;
	5 6 7
]

# ╔═╡ 7ac87c2a-52a8-4e5f-be5d-517d68adc6af
big_determinant(q1B) # Not independent

# ╔═╡ d1d81554-a38e-422c-a0ac-762e98c3f9a3
cofactor_determinant(q1B)

# ╔═╡ a09e4c71-a672-41f8-8989-14f583da0006
q1C = [
	1 1 1;
	1 1 0;
	1 0 0
]

# ╔═╡ bfbc800c-70ca-4bc0-8e85-db7b0c5d0ca2
big_determinant(q1C)

# ╔═╡ a8aff092-254d-416d-a6f3-1e7a6ddbe1a5
cofactor_determinant(q1C)

# ╔═╡ de645e81-575a-416a-ade9-728ac5f915bd
md"""
## Q2
"""

# ╔═╡ 60ada824-3740-40c0-bc55-0946283db92e
q2A = [
	1 1 0;
	1 0 1;
	0 1 1
]

# ╔═╡ f023ce25-7236-482c-bb89-6e77f9444c60
big_determinant(q2A)

# ╔═╡ 660eea34-7f99-4214-bccd-1ee0a299405a
q2B = [
	1 2 3;
	4 5 6;
	7 8 9
]

# ╔═╡ decdbe4a-3239-4930-b014-8d645477e6af
big_determinant(q2B) # Not independent

# ╔═╡ 61948f65-5a96-415d-88f9-d05169b7d775
q2C = Int.(vcat(hcat(q2A,zeros(3, 3)), hcat(zeros(3,3),q2A)))

# ╔═╡ 3a64634c-1afa-42f5-bed5-419b3eaf1a52
big_determinant(q2C)

# ╔═╡ c73c3aa1-f016-48e7-9dae-2ca9d2e71369
det(q2C)

# ╔═╡ 478ce283-3b5c-4a19-8aac-a8e7678928b0
q2D = Int.(vcat(hcat(q2A,zeros(3, 3)), hcat(zeros(3,3),q2B)))

# ╔═╡ 3a3f3a0d-8ecf-492e-a0a7-34f87cf35842
big_determinant(q2D)

# ╔═╡ 6a5a0108-6fec-413f-a063-6e301d621f31
det(q2D)

# ╔═╡ c4b032d4-9ab4-4ebd-baed-8fe2ca0a9a64
md"""
## Q3
"""

# ╔═╡ a041250a-6a9a-47d8-9dbe-ee562722a53f
begin
	q3A = [
		x x x;
		0 0 x;
		0 0 x
	]
	show_num(q3A)
end

# ╔═╡ b357c956-214e-42e1-9af9-d918b34498e1
show_num(big_determinant(q3A))

# ╔═╡ 5c6392f6-f4c3-40d0-b029-2336c34e2e37
show_num(cofactors(q3A))

# ╔═╡ 24c32878-79d1-4747-87db-b79e7d39b7a6
show_num(lu(q3A, check=false).U)  # Rank < 3

# ╔═╡ 3629a2a9-246e-48f0-86e7-433d2fc28a29
show_num([q3A[1,j]*(((-1)^(1 + j)) * cofactor_determinant(cofactor_matrix(q3A,1,j))) for j=1:size(q3A)[2]])

# ╔═╡ e2d5113f-16a5-4e6e-ad30-8965f78192d2
show_num(big_determinant_terms(q3A))

# ╔═╡ daf13d87-4ff6-4cdb-8422-c43eb426ac57
md"""
## Q4
"""

# ╔═╡ 263e7145-0eb5-4e43-a5d5-5a714791ef1d
q4A = [
	1 0 0 1;
	0 1 1 1;
	1 1 0 1;
	1 0 0 1
]

# ╔═╡ dd25c48b-182a-4ad2-8da9-18c73950b2c7
q4P1 = [
	1 0 0 0;
	0 0 1 0;
	0 1 0 0;
	0 0 0 1
]

# ╔═╡ 41b92699-0134-451c-a1c2-e44b4132e41b
q4P2 = [
	0 0 0 1;
	0 0 1 0;
	0 1 0 0;
	1 0 0 0
]

# ╔═╡ 3d512797-1045-4c78-9a17-c6aa0d7a80a9
q4B = [
	1 0 0 2;
	0 3 4 5;
	5 4 0 3;
	2 0 0 1
]

# ╔═╡ 8850715b-6e5d-474c-9623-0c498e1722c8
big_determinant(q4A) # 1 -1

# ╔═╡ 788798b7-775c-428b-a3ca-78ca03a10d04
big_determinant(q4B)

# ╔═╡ 672ce5a4-0f6a-48ce-bb54-dad1371ea60c
md"""
## Q5
"""

# ╔═╡ 6cc55784-471a-4c66-a648-243f0707d6a5

q5A = [
	a b c d;
	e f g h;
	i j k l;
	0 0 0 0
];

# ╔═╡ 38fc0774-6ffe-4ca1-874a-e9b3fcf79b82
show_num(big_determinant(q5A))

# ╔═╡ 0e17177d-4e24-49f1-8a29-2620f7108362
md"""
I(4) has the most possible 0s 
"""

# ╔═╡ 03b9475f-efdf-4897-8797-e8498f343062
md"""
## Q6

a)
``\det(A) = \sum{\det(P)a_{1\alpha}a_{2\beta}a_{3\gamma}}; \alpha \neq 1, \beta \neq 2, \gamma \neq 3``

``\det(A) = a_{12}a_{23}a_{32} + a_{13}a_{21}a_{32}``

"""

# ╔═╡ 9f7feabb-48b1-4285-9fca-4ea539036d91
q6Aa = [
	0 a12 a13;
	a21 0 a23;
	a31 a32 0;
];

# ╔═╡ 0e56b641-6332-4c2f-af47-e4bef2fdf87a
show_num(big_determinant_terms(q6Aa))

# ╔═╡ c851188c-fe5e-4cd9-a281-d51fa3ee6231
q6Ab = [
	0 a12 a13 a14;
	a21 0 a23 a24;
	a31 a32 0 a34;
	a41 a42 a43 0;
];

# ╔═╡ 14a3fdb9-af9e-40cf-bc07-43838fc352b1
show_num(big_determinant_terms(q6Ab))

# ╔═╡ d2f3b244-dc7f-4a51-99bb-cfeffb64d21b
md"""
## Q7
"""

# ╔═╡ 44218092-fd43-4ed2-9861-77e9d77e8500
q7P = [Permutation(5, i) for i=1:factorial(5) if det(Matrix(Permutation(5, i))) == 1]

# ╔═╡ 9a5416e8-4dc7-48ba-a2fb-71c6fcae914f
length(q7P)

# ╔═╡ 835fa8bd-3c62-4e43-9792-e5d8a5ef3d96
Matrix(q7P[end])

# ╔═╡ 9361ae3c-ab6a-4c9b-8eee-bac81c78eb98
md"""
## Q8

If there is a non-zero determinant term, you can simply use the permutation operation associated with the term to get values of hte diagonal. 
"""

# ╔═╡ c3bfaea3-db75-456a-ae17-9bdfc5cacdef
md"""
## Q9

Hadamard's maximum determinant game!
"""

# ╔═╡ 6dce6c59-215e-4838-8263-e0bcdfce5f4b
function smart_maximal_determinant(q33size)
	reduced_size = q33size - 1
	max_determinant = 0
	max_determinant_matrix = zeros(reduced_size, reduced_size)
	show(max_determinant)

	for q33_value = 1:2^reduced_size^2
		q33_matrix = reshape(map(x-> x=='0' ? 0 : 1, collect(bitstring(UInt64(q33_value))[64-(reduced_size^2-1):end])), (reduced_size,reduced_size))
		if det(q33_matrix) > max_determinant
			max_determinant = det(q33_matrix)
			max_determinant_matrix = q33_matrix
		end
	end
	
	max_determinant * 2^reduced_size, vcat(ones(1, q33size), hcat(ones(reduced_size, 1), map(x-> x == 0 ? 1 : -1, max_determinant_matrix))) 
end

# ╔═╡ 9ed209f0-cb59-44c3-bcd2-17a0d984204c
smart_maximal_determinant(3)

# ╔═╡ 6a133cee-048e-4c26-9c11-be71a22edddc
md"""
## Q10
"""

# ╔═╡ 27900e9d-3ebe-4cc8-87f7-96573688f6ca
q10Ps = [Matrix(Permutation(4, i)) for i=1:factorial(4) if det(Matrix(Permutation(4, i))) == 1]

# ╔═╡ ae9b324d-b066-4b43-b165-170939c475e3
size(q10Ps)

# ╔═╡ 7cf14733-af0e-4546-8f54-1dc51662e3cf
det.(map(x->x +I, q10Ps))

# ╔═╡ fbe30313-f9e9-4151-b426-df61b56dc28f
md"""
## Q11
"""

# ╔═╡ 433b3e27-fd82-4b19-ad90-102618fcc381
q11A = [
	a b;
	c d
];

# ╔═╡ aaec4cc7-817d-41fa-9b50-3f9f90ffeb27
q11C = all_cofactors(q11A);

# ╔═╡ 7627d69d-2ae2-420f-96c1-64f740aef884
show_num(q11C)

# ╔═╡ 939e461c-82fd-4136-ab07-904facfad3f4
show_num(q11A*vcat(q11C))

# ╔═╡ eddef1e3-ffff-400e-af18-88cfe67a81b2
q11B = [
	1 2 3;
	4 5 6;
	7 0 0
]

# ╔═╡ dac6fdd4-315e-40ea-8949-14e2e7bcac84
q11D = all_cofactors(q11B)

# ╔═╡ 30c4f07b-3d14-4fcc-99a9-e6aca92da1cf
cofactor_determinant(q11B)

# ╔═╡ 339633aa-1eb7-42b4-b922-6a57fcc57455
md"""
## Q12
"""

# ╔═╡ 68a2d75c-04da-4a24-b823-d1ec781289c1
q12A = [
	2 -1 0;
	-1 2 -1;
	0 -1 2
]

# ╔═╡ 269120fd-b9b2-4e90-9c67-c83a0dc7d779
rat.(inv(q12A))

# ╔═╡ 6226ca5f-0ed3-4e63-86c0-79f561db197e
q12CT = all_cofactors(q12A)

# ╔═╡ 57d77d29-9048-4f95-9d7c-5d521dd3bd81
q12A * q12CT'

# ╔═╡ 21f0358e-bf22-4153-9c0b-eed34a51722b
rat.((1/det(q12A))*q12CT')

# ╔═╡ 9a7b4489-1f5a-4589-9d93-d26fcc470dca
md"""
## Q13
"""

# ╔═╡ 23dd8a67-dbf5-4870-b278-d34336dfac30
q13C1 = [0]

# ╔═╡ 8f77d447-23c2-4832-95f1-c66919e1f3d7
q13C2 = [
	0 1;
	1 0
]

# ╔═╡ 208af161-d36c-43b7-adc8-16ad987b8de7
tridiagonal13(size) = Tridiagonal(ones(size-1), zeros(size), ones(size-1))

# ╔═╡ 26ad19ea-6d75-4ae4-9aba-269312120963
tridiagonal13(2)

# ╔═╡ 32aa3e55-042e-4ab1-8273-7291a8de01e6
det(tridiagonal13(1))

# ╔═╡ 58bd585b-ac46-40d2-b51f-10c0f56225d2
q13C3 = tridiagonal13(3)

# ╔═╡ e4570afd-c1b1-49e5-a4a8-29a092a45c61
q13C4 = tridiagonal13(4)

# ╔═╡ 2119b87a-4293-4e89-af68-d64d3afe7215
[big_determinant(Array(tridiagonal13(i))) for i=1:10]

# ╔═╡ 1c7be823-cae7-4f7f-bcf2-7c56af951e38
md"""
The bottom corner of the right corner of ``C_n = -C_{n-2}``
"""

# ╔═╡ b4dde61d-44fc-42ed-9e33-4efcbfcc9bc3
md"""
## Q14

The relationship is ``sign(C_n) = (-1)^{\frac{n}{2}}``
"""

# ╔═╡ 221f3f2d-442e-48bd-aa2c-5fba5e069844
md"""
## Q15


Notice that the cofactor equation is ``\det(A) = a_{11}E_{11} + a_{12}E_{12} + ...``

Also notice that ``E_{12}`` while be negative because ``(-1)^{(i+j)}``

But because ``E`` is ``(1,1,1)``-tridiagonal, only the first two cofactors are needed, so: 

``E_n = E_{n-1} - E_{n-2}``
"""

# ╔═╡ 44fff8ac-9940-4cb4-9ab2-ee5ea9ced918
tridiagonal15(size) = Tridiagonal(ones(size-1), ones(size), ones(size-1))

# ╔═╡ fa0908f6-56d6-4139-9978-c6621d3d8430
[cofactor_determinant(tridiagonal15(i)) for i=1:8]

# ╔═╡ 75b5757f-b3f4-41bb-90a6-feb50d99e2d3
@memoize function fast_15_determinant(i)
	if i == 1
		return 1
	elseif i == 2
		return 0
	else
		return fast_15_determinant(i-1) - fast_15_determinant(i-2)
	end
end

# ╔═╡ 6d621d9c-bea6-4b74-9f8c-3851a4026e3a
[fast_15_determinant(i) for i=1:8]

# ╔═╡ 7f318510-c5f2-4f69-b638-62e93e270dfb
fast_15_determinant(100)

# ╔═╡ c5a41df2-b049-4ec7-947d-22c97ef65b38
md"""
## Q16

By the same argument as (15), only the first two cofactors are needed to specify the deteriminant through cofactors: ``\det(A) = a_{11}E_{11} + a_{12}E_{12}``, also the negativity of ``E_{12}`` cancels out the permutation negativity.
"""

# ╔═╡ 5eeb5374-7097-4dce-9921-23871382b049
tridiagonal16(size) = Tridiagonal(ones(size-1), ones(size), -1*ones(size-1))

# ╔═╡ 881f05dc-f9c8-4b91-95d8-abd429c06aff
q16F2 = tridiagonal16(2)

# ╔═╡ 92b5db4d-0194-4a96-80a7-e6bbea11e5c1
q16F3 = tridiagonal16(3)

# ╔═╡ 8dceb14a-4513-42f8-b502-063b158ebe41
q16F4 = tridiagonal16(4)

# ╔═╡ f939ba5d-7f82-48db-a830-4b4c487b519d
md"""
## Q17

``|B_n| = 2|B_{n-1}``
"""

# ╔═╡ 292f3d65-7669-4326-9597-73313345c3e9
tridiagonal17(size) = Tridiagonal(-1*ones(size-1), vec(hcat(1, ones(size-1)*2...)), -1*ones(size-1))

# ╔═╡ 25e06f36-f583-41c1-8aca-a92bcc780dbf
q17B4 = tridiagonal17(4)

# ╔═╡ 64d53601-3b03-4155-b609-b3864fcc5b31
q17B3 = tridiagonal17(3)

# ╔═╡ 14ecaf00-2497-45b6-9c38-65644ff59745
q17B2 = tridiagonal17(2)

# ╔═╡ e22425ea-35fa-4a70-a7db-b081d9062c80
cofactor_determinant(cofactor_matrix(q17B4, 4, 4))

# ╔═╡ 51c925ed-5df6-427e-b790-871c2f198bb9
[cofactor_determinant(tridiagonal17(i)) for i=2:4]

# ╔═╡ 07105032-916f-4caa-850f-fdacd8a1b971
md"""
## Q18
"""

# ╔═╡ 21185fd0-2d71-4293-82d0-e86e63e6c9f8
tridiagonal18a(size) = Tridiagonal(-1*ones(size-1), 2*ones(size), -1*ones(size-1))

# ╔═╡ ace0dd4d-d856-4468-a6c6-2b5e23c4f991
tridiagonal18b(size) = Tridiagonal(-1*ones(size-1), vec(hcat(1, ones(size-1)*2...)), vec(hcat(0.0, -1*ones(size-2)...)))

# ╔═╡ 3c8eeb2a-8b60-4922-bbf3-ae2910b2870a
tridiagonal18b(2)

# ╔═╡ 2aae9bed-9fa2-467e-86bc-284acb8defe4
[cofactor_determinant(tridiagonal18a(i)) - cofactor_determinant(tridiagonal18a(i - 1)) for i=2:5]

# ╔═╡ c3428cf4-cf48-4b9d-ba3c-da8187018bdd
md"""
## Q19
"""

# ╔═╡ e2c7912a-d980-466a-8843-ae13893839b6
q19V4 = [
	1 a a^2 a^3;
	1 b b^2 b^3;
	1 c c^2 c^3;
	1 x x^2 x^3
];

# ╔═╡ 185c88a1-5f39-4459-a85a-3e91c025d242
show_num(simplify(cofactor_determinant(substitute.(q19V4, (Dict([x => a]),))), expand=true))

# ╔═╡ 86e8d8bf-6f0d-4f6f-b46d-1cb648a93af3
show_num(simplify(cofactor_determinant(substitute.(q19V4, (Dict([x => b]),))), expand=true))

# ╔═╡ 9562bb5f-cd86-4e6c-a9f7-bff4c5014067
show_num(simplify(cofactor_determinant(substitute.(q19V4, (Dict([x => c]),))), expand=true))

# ╔═╡ 617dbdac-1bfa-4834-9d7a-19340204f5c0


# ╔═╡ 6d0194a4-ff29-49dd-8146-1de35f77998b
md"""
## Q20
"""

# ╔═╡ fbe51c6c-5d45-4204-b7e3-d43bb09ecc8b
constructq20(size) = ones(size, size) - I

# ╔═╡ 6ae1dd18-b181-4d9d-8246-167cf9535df0
q20G2 = constructq20(2)

# ╔═╡ 7f261229-8220-4be5-980e-0141c3cc31e1
q20G3 = constructq20(3)

# ╔═╡ 86fd8e42-2abc-4b67-9950-1480f7ba636d
q20G4 = constructq20(4)

# ╔═╡ e7cffa17-9a2b-4093-9f01-08f4443ecef3
[det(constructq20(i)) for i=1:4]

# ╔═╡ 88a6af52-c9b6-478a-8b80-964e9d54fc55
md"""
I predict that ``G_n = (-1)^{(n-1)}*(n-1)``  
"""

# ╔═╡ e7837ea4-f805-45e1-9534-347d253a2de3
md"""
## Q21
"""

# ╔═╡ 50aee535-4ec6-467f-ac00-9902b3ead43c
constructq21(size) = Tridiagonal(ones(size-1), 3*ones(size), ones(size-1))

# ╔═╡ 99e8e7a2-8591-4f0c-8830-07cc0883b0b6
q21S1 = [3]

# ╔═╡ b9172bca-f008-4579-b4c4-d1b7fc1e9a41
q21S2 = constructq21(2)

# ╔═╡ 3235d955-ad52-4afd-8eb4-e7f6eb9dd03f
q21S3 = constructq21(3)

# ╔═╡ 80d43e49-eee3-42b6-ba62-76c9671b9f65
[det(constructq21(i)) for i=1:4]

# ╔═╡ 84449c10-438a-44b8-ac1e-dd617296dc8a
md"""
## Q22
"""

# ╔═╡ e4eb0ee4-76eb-4b63-ae81-6efaf90c5338
constructq22(size) = Tridiagonal(ones(size-1), vec(hcat(2.0, ones(size-1)*3...)), ones(size-1))

# ╔═╡ 36dd5725-1589-4f29-ae6c-9f4ebe223e83
[cofactor_determinant(constructq22(i)) for i=1:5]

# ╔═╡ 7379277e-aec8-4bc6-b763-cadd0ebd6314
md"""
Recall that the cofactor determinant formula is:

``\det(A) = \sum_{i=1}^n a_{1i}C_{1i}``

For this matrix the equation is:

``\det(A) = 2|M_{11}| - |M_{12}|``

``|M_{11}| = |S_{n-1}|``

``|M_{12}| = |S_{n-2}|``

By replacing the ``A[1, 1]` = 3`` with a ``2``, the determinant equation went from.

``\det(A) = 3|S_{n-1}| - |S_{n-2}|`` to ``\det(A) = 2|S_{n-1}| - |S_{n-2}|``

"""

# ╔═╡ a4b5f824-e256-40f8-8d59-02605f0954f7
[det(cofactor_matrix(constructq22(i), 1, 2)) for i=2:6]

# ╔═╡ 4c528352-37d6-4cdb-8837-b4234b526f68
[det(cofactor_matrix(constructq22(i), 1, 1)) for i=2:6]

# ╔═╡ b9ebfa65-f18a-47c5-a4e6-4af7305d9261
constructq22(3)

# ╔═╡ a886d29e-bce2-4316-860d-3063df707f36
md"""
## Q23

a) The statement is true because of the formula for the 2x2 deteriment:

``|U| = |A||D| - |B||C|``, but since ``|C|=0``:
``|U| = |A||D|``


"""

# ╔═╡ fabb97a8-701c-4d37-b679-a254b95b0a26
q23A = [
	1 2;
	3 4
]

# ╔═╡ 015073d4-9e93-40df-b69b-3fd7b0026c72
q23B = [
	1 -1;
	2 7
]

# ╔═╡ ac0da309-8776-45b6-84e3-cdd8f574ac7c
q23C = [
	6 10;
	7 4
]

# ╔═╡ c587da33-bd3a-4ed2-8cdc-61cbcb7b86d9
q23D = [
	5 3;
	8 4
]

# ╔═╡ 075e3263-a63f-4172-aacb-ed47ef4b1e79
rat.(det(vcat(hcat(q23A, q23B), hcat(q23C, q23D))))

# ╔═╡ 101dcdc7-b40b-4f29-aafd-451ecfea8ac7
rat.(det(q23A)*det(q23D) - det(q23C)*det(q23B))

# ╔═╡ 01059502-8b67-41be-9a89-fd2fbb938267
rat.(det(q23A*q23D - q23C*q23B))

# ╔═╡ 4447861d-56f6-44ae-91d6-87a393be0a43
md"""
## Q24

a) All ``L``- pivots will be ``1``.

General formula:

``d_k*det(A_{k-1}) = det(A_k)``

``|A_1| = |U_1| = 2``

``|A_2| = |U_2| = 3*2 = 6``

``|A_3| = |U_3| =  -1*6 = -6``

b) 
``|A_1| = 5``
``|A_2| = \frac{6}{5}``
``|A_3| = \frac{7}{6}``

"""

# ╔═╡ eb0c5ef2-28d5-431f-8d22-7db0564c24a3
md"""
## Q25

``\v{I & 0 \\ -CA^{-1} & I} = I^2 - 0 = I``

``\v{A & B \\ 0 & D-CA^{-1}B} = |A||D - CA^{-1}B| = |AD -ACA^{-1}B| = |AD-CB|``
"""


# ╔═╡ 6a6ec27c-4cd1-429b-b881-f92552550118
md"""
## Q26
"""

# ╔═╡ 76566bd0-00b0-4ce4-871a-ed6478dd03f1
q26M1 = [
	0 a1 a2;
	b1 1 0
	b2 0 1
];

# ╔═╡ 3a70587b-4902-471c-8516-5fa79d2c3e4f
show_num(big_determinant(q26M1))

# ╔═╡ 2278063b-9608-4215-8d7e-fbeb948ae6dc
q26M2 = [
	0 0 a1;
	0 0 a2;
	b1 b2 1
];

# ╔═╡ 0ad4b081-015c-4544-8db0-1e1f2e6f8c22
show_num(big_determinant(q26M2))

# ╔═╡ 63e2463d-fa27-44cf-a2cc-fc53caa019f6
md"""
## Q27

``\frac{\partial}{\partial a_{11}} \det(A) = \frac{\partial}{\partial a_{11}} a_{11}C_{11} + ... = C_{11}``
"""

# ╔═╡ 13b14998-2bf9-4e5d-b902-d480a5d9711c
md"""
## Q28
"""

# ╔═╡ f1cda707-f1c6-479c-80a7-fdf13e818528
q28D = [
	1 2 3;
	4 5 6;
	7 8 9
]

# ╔═╡ 9168fe26-f50d-42c5-9947-8b7af9701257
q28diags = [sign(Permutation(3, i))*reduce(*, diag(Matrix(Permutation(3, i))*q28D)) for i=1:6]

# ╔═╡ dd6ebd84-d654-4b85-94cf-5f08be08e2ac
reduce(+,q28diags) 

# ╔═╡ 2b9ea586-69c2-41b8-97f3-2be58186652b
md"""
These values all add up to 0, which suggests that this matrix is singular. Also ``2R_2 - R_1 = R_3``
"""

# ╔═╡ 24ddccd5-8787-4c58-8c62-4dca6787762f
md"""
## Q29
"""

# ╔═╡ 83476f5b-9e21-40e0-a91f-32f07e42af9a
q29_terms = big_determinant_terms(Array(tridiagonal15(4)))

# ╔═╡ 422ba08d-fea2-41fe-b599-0ba01fe1f0a6
sum(q29_terms)

# ╔═╡ 03399151-21e9-4e5f-aa6e-a7813ca47a15
md"""
## Q30
"""

# ╔═╡ cb514a8a-123c-4839-a17d-8d1fb41430ab
constructq30(size) = Tridiagonal(ones(size-1)*-1, ones(size) * 2, ones(size-1) * -1)

# ╔═╡ dd6451df-c02e-4511-ab39-fcf5c1075cd2
[Permutation(4, i) for i=1:factorial(4) if big_determinant_terms(Array(constructq30(4)))[i] != 0]

# ╔═╡ 59d66bf5-20d9-4eba-b245-28341aec58c0
md"""
## P31
"""

# ╔═╡ 0be7d6ab-313c-4251-a3fc-014d0ff08d8e
q31P = [
	0 0 0 1;
	1 0 0 0;
	0 1 0 0;
	0 0 1 0
]

# ╔═╡ a83b79a2-176a-4055-9e85-c7ac48b80426
cofactor_determinant(q31P)

# ╔═╡ fa662437-2904-4eaa-b686-ad690d5775e0
big_determinant(q31P)

# ╔═╡ 50d8b353-12d0-45f2-833a-0199f18b7b1f
md"""
3 exchanges are required to reorder ``4, 1, 2, 3``
"""

# ╔═╡ 90f6ddd0-76b8-442e-a976-84a0f8233c7c
cofactor_determinant(q31P^2)

# ╔═╡ 8baed9ad-8824-41f8-ab8a-c2acf60edb46
md"""
## Q32

``k = 2n +2``

``F_k = F_{k-1} + F_{k-2}``

``F_{2n+2} = F_{2n+1} + F_{2n}``

``F_{2n+1} = F_{2n} + F_{2n-1}``

``F_{2n+2} = 2F_{2n} + F_{2n-1}``

By rearranging the Fibonacci identity:

``F_{k} - F_{k-2} = F_{k-1}``

If ``k = 2n-1``

``F_{2n+2} = 3F_{2n} - F_{2n-2} \blacksquare``

"""

# ╔═╡ 0a8fad3e-759e-4b8b-a3a0-466574feb9c9
md"""
## Q33
"""

# ╔═╡ 4f527a3c-5e0d-4b7d-916f-ff580eb8749e
q33A = matrixdepot("pascal", 4)

# ╔═╡ 1cd02f96-6e4e-44d9-b44d-a15147b626ba
begin
	q33Z = zeros(4,4)
	q33Z[4,4] = 1
end

# ╔═╡ fad5b501-185b-4a0a-863c-2fcebb933385
cofactors(q33A-q33Z)

# ╔═╡ 6df45b5a-69b1-45d9-985d-72cd6cd4f55e
sum(cofactors(q33A-q33Z)) # The cofactors add to 0!

# ╔═╡ f8653bd3-d2f0-43fd-9d64-381aa7a5c5aa
md"""
## Q34
"""

# ╔═╡ dfe40c17-4f5d-4a63-b6d1-cb4d6f4d178d
q34A = [
	x x x x x;
	x x x x x;
	0 0 0 x x;
	0 0 0 x x;
	0 0 0 x x
];

# ╔═╡ 9254e969-0f77-4cdb-9198-618fee92d1db
md"""
## Q35
"""

# ╔═╡ Cell order:
# ╠═8993d71d-9c83-4b2b-87be-fe5afc2fd359
# ╠═45411a2e-449a-4475-838b-187adabfc6b1
# ╠═60bb7f84-877a-4e32-b50e-8c63a5c2b24b
# ╠═16717b6f-4e60-41d6-b796-f2b1c09decc8
# ╠═aa4a53fa-7168-408b-92ed-6c6bb4a4d128
# ╠═48f45c84-2fae-4c89-9d1d-b5fdff136ffd
# ╠═718e651f-a9e6-4e6d-ab50-9a429179bdf8
# ╠═3b1025b8-ab31-41d2-9ff5-1d31e7e1ecf2
# ╠═111fa6a0-3d84-420e-b546-59075a03a5c2
# ╠═e4c74c4c-97c4-4865-bc54-cddf00fa2e58
# ╠═c3e4419b-7421-423f-a804-0b98f7638598
# ╠═742113f5-7cca-4139-b288-a0dc41264423
# ╠═5b4a23b6-17ee-4ce9-814f-cb8dd2489a58
# ╠═49f2380f-af3b-4b20-afde-1b793499cdec
# ╠═e9ed900c-a286-4e29-adaf-5fff401068f9
# ╠═f880b525-d950-472a-a966-9313c3ad9e14
# ╠═95074837-fbab-43ba-812c-680895172946
# ╠═e210ddda-0ebc-42fe-88b5-6545838d54e1
# ╠═fa7f4e03-f93e-401a-9280-0f94957d22ac
# ╠═400ee06e-1172-4805-abf9-981a2a66dadf
# ╠═b93f597c-fae2-4ed2-ae24-b2f17a94d8c8
# ╠═7ae07d75-d7e4-467c-8a0a-da084de0d5ac
# ╠═aa4721fa-046c-4d28-a8f6-245d4e73a7d0
# ╠═862e32ee-e26a-4cd5-ae0f-83a24b898102
# ╠═1b83dbef-8bc2-4fe0-a357-7f95b9a2d455
# ╠═9d4a50b6-8fe5-4611-aafe-af548bc5970d
# ╠═95050319-88ab-4c0c-9406-d1dcc36c2ccd
# ╠═2b5d4b30-aeab-11eb-0b6d-5521705ce8f5
# ╠═b18b009f-f2ea-4093-a91c-1bf56c661a86
# ╠═68fe8dde-2fb3-46b6-9d4f-cc9f85c78ee2
# ╠═642bef99-6b5f-41c3-84c0-98a47ab53f21
# ╠═6cce2d8f-2d6d-46e4-917b-a06078acc12a
# ╠═33714c08-90c6-4044-ad80-c5649024f47f
# ╠═7ac87c2a-52a8-4e5f-be5d-517d68adc6af
# ╠═d1d81554-a38e-422c-a0ac-762e98c3f9a3
# ╠═a09e4c71-a672-41f8-8989-14f583da0006
# ╠═bfbc800c-70ca-4bc0-8e85-db7b0c5d0ca2
# ╠═a8aff092-254d-416d-a6f3-1e7a6ddbe1a5
# ╠═de645e81-575a-416a-ade9-728ac5f915bd
# ╠═60ada824-3740-40c0-bc55-0946283db92e
# ╠═f023ce25-7236-482c-bb89-6e77f9444c60
# ╠═660eea34-7f99-4214-bccd-1ee0a299405a
# ╠═decdbe4a-3239-4930-b014-8d645477e6af
# ╠═61948f65-5a96-415d-88f9-d05169b7d775
# ╠═3a64634c-1afa-42f5-bed5-419b3eaf1a52
# ╠═c73c3aa1-f016-48e7-9dae-2ca9d2e71369
# ╠═478ce283-3b5c-4a19-8aac-a8e7678928b0
# ╠═3a3f3a0d-8ecf-492e-a0a7-34f87cf35842
# ╠═6a5a0108-6fec-413f-a063-6e301d621f31
# ╠═c4b032d4-9ab4-4ebd-baed-8fe2ca0a9a64
# ╠═a041250a-6a9a-47d8-9dbe-ee562722a53f
# ╠═b357c956-214e-42e1-9af9-d918b34498e1
# ╠═5c6392f6-f4c3-40d0-b029-2336c34e2e37
# ╠═24c32878-79d1-4747-87db-b79e7d39b7a6
# ╠═3629a2a9-246e-48f0-86e7-433d2fc28a29
# ╠═e2d5113f-16a5-4e6e-ad30-8965f78192d2
# ╠═daf13d87-4ff6-4cdb-8422-c43eb426ac57
# ╠═263e7145-0eb5-4e43-a5d5-5a714791ef1d
# ╠═dd25c48b-182a-4ad2-8da9-18c73950b2c7
# ╠═41b92699-0134-451c-a1c2-e44b4132e41b
# ╠═3d512797-1045-4c78-9a17-c6aa0d7a80a9
# ╠═8850715b-6e5d-474c-9623-0c498e1722c8
# ╠═788798b7-775c-428b-a3ca-78ca03a10d04
# ╠═672ce5a4-0f6a-48ce-bb54-dad1371ea60c
# ╠═6cc55784-471a-4c66-a648-243f0707d6a5
# ╠═38fc0774-6ffe-4ca1-874a-e9b3fcf79b82
# ╠═0e17177d-4e24-49f1-8a29-2620f7108362
# ╠═03b9475f-efdf-4897-8797-e8498f343062
# ╠═9f7feabb-48b1-4285-9fca-4ea539036d91
# ╠═0e56b641-6332-4c2f-af47-e4bef2fdf87a
# ╠═c851188c-fe5e-4cd9-a281-d51fa3ee6231
# ╠═14a3fdb9-af9e-40cf-bc07-43838fc352b1
# ╠═d2f3b244-dc7f-4a51-99bb-cfeffb64d21b
# ╠═44218092-fd43-4ed2-9861-77e9d77e8500
# ╠═9a5416e8-4dc7-48ba-a2fb-71c6fcae914f
# ╠═835fa8bd-3c62-4e43-9792-e5d8a5ef3d96
# ╠═9361ae3c-ab6a-4c9b-8eee-bac81c78eb98
# ╠═c3bfaea3-db75-456a-ae17-9bdfc5cacdef
# ╠═6dce6c59-215e-4838-8263-e0bcdfce5f4b
# ╠═9ed209f0-cb59-44c3-bcd2-17a0d984204c
# ╠═6a133cee-048e-4c26-9c11-be71a22edddc
# ╠═27900e9d-3ebe-4cc8-87f7-96573688f6ca
# ╠═ae9b324d-b066-4b43-b165-170939c475e3
# ╠═7cf14733-af0e-4546-8f54-1dc51662e3cf
# ╠═fbe30313-f9e9-4151-b426-df61b56dc28f
# ╠═433b3e27-fd82-4b19-ad90-102618fcc381
# ╠═aaec4cc7-817d-41fa-9b50-3f9f90ffeb27
# ╠═7627d69d-2ae2-420f-96c1-64f740aef884
# ╠═939e461c-82fd-4136-ab07-904facfad3f4
# ╠═eddef1e3-ffff-400e-af18-88cfe67a81b2
# ╠═dac6fdd4-315e-40ea-8949-14e2e7bcac84
# ╠═30c4f07b-3d14-4fcc-99a9-e6aca92da1cf
# ╠═339633aa-1eb7-42b4-b922-6a57fcc57455
# ╠═68a2d75c-04da-4a24-b823-d1ec781289c1
# ╠═269120fd-b9b2-4e90-9c67-c83a0dc7d779
# ╠═6226ca5f-0ed3-4e63-86c0-79f561db197e
# ╠═57d77d29-9048-4f95-9d7c-5d521dd3bd81
# ╠═21f0358e-bf22-4153-9c0b-eed34a51722b
# ╠═9a7b4489-1f5a-4589-9d93-d26fcc470dca
# ╠═23dd8a67-dbf5-4870-b278-d34336dfac30
# ╠═8f77d447-23c2-4832-95f1-c66919e1f3d7
# ╠═208af161-d36c-43b7-adc8-16ad987b8de7
# ╠═26ad19ea-6d75-4ae4-9aba-269312120963
# ╠═32aa3e55-042e-4ab1-8273-7291a8de01e6
# ╠═58bd585b-ac46-40d2-b51f-10c0f56225d2
# ╠═e4570afd-c1b1-49e5-a4a8-29a092a45c61
# ╠═2119b87a-4293-4e89-af68-d64d3afe7215
# ╠═1c7be823-cae7-4f7f-bcf2-7c56af951e38
# ╠═b4dde61d-44fc-42ed-9e33-4efcbfcc9bc3
# ╠═221f3f2d-442e-48bd-aa2c-5fba5e069844
# ╠═44fff8ac-9940-4cb4-9ab2-ee5ea9ced918
# ╠═fa0908f6-56d6-4139-9978-c6621d3d8430
# ╠═75b5757f-b3f4-41bb-90a6-feb50d99e2d3
# ╠═6d621d9c-bea6-4b74-9f8c-3851a4026e3a
# ╠═7f318510-c5f2-4f69-b638-62e93e270dfb
# ╠═c5a41df2-b049-4ec7-947d-22c97ef65b38
# ╠═5eeb5374-7097-4dce-9921-23871382b049
# ╠═881f05dc-f9c8-4b91-95d8-abd429c06aff
# ╠═92b5db4d-0194-4a96-80a7-e6bbea11e5c1
# ╠═8dceb14a-4513-42f8-b502-063b158ebe41
# ╠═f939ba5d-7f82-48db-a830-4b4c487b519d
# ╠═292f3d65-7669-4326-9597-73313345c3e9
# ╠═25e06f36-f583-41c1-8aca-a92bcc780dbf
# ╠═64d53601-3b03-4155-b609-b3864fcc5b31
# ╠═14ecaf00-2497-45b6-9c38-65644ff59745
# ╠═e22425ea-35fa-4a70-a7db-b081d9062c80
# ╠═51c925ed-5df6-427e-b790-871c2f198bb9
# ╠═07105032-916f-4caa-850f-fdacd8a1b971
# ╠═21185fd0-2d71-4293-82d0-e86e63e6c9f8
# ╠═ace0dd4d-d856-4468-a6c6-2b5e23c4f991
# ╠═3c8eeb2a-8b60-4922-bbf3-ae2910b2870a
# ╠═2aae9bed-9fa2-467e-86bc-284acb8defe4
# ╠═c3428cf4-cf48-4b9d-ba3c-da8187018bdd
# ╠═e2c7912a-d980-466a-8843-ae13893839b6
# ╠═185c88a1-5f39-4459-a85a-3e91c025d242
# ╠═86e8d8bf-6f0d-4f6f-b46d-1cb648a93af3
# ╠═9562bb5f-cd86-4e6c-a9f7-bff4c5014067
# ╠═617dbdac-1bfa-4834-9d7a-19340204f5c0
# ╠═6d0194a4-ff29-49dd-8146-1de35f77998b
# ╠═fbe51c6c-5d45-4204-b7e3-d43bb09ecc8b
# ╠═6ae1dd18-b181-4d9d-8246-167cf9535df0
# ╠═7f261229-8220-4be5-980e-0141c3cc31e1
# ╠═86fd8e42-2abc-4b67-9950-1480f7ba636d
# ╠═e7cffa17-9a2b-4093-9f01-08f4443ecef3
# ╠═88a6af52-c9b6-478a-8b80-964e9d54fc55
# ╠═e7837ea4-f805-45e1-9534-347d253a2de3
# ╠═50aee535-4ec6-467f-ac00-9902b3ead43c
# ╠═99e8e7a2-8591-4f0c-8830-07cc0883b0b6
# ╠═b9172bca-f008-4579-b4c4-d1b7fc1e9a41
# ╠═3235d955-ad52-4afd-8eb4-e7f6eb9dd03f
# ╠═80d43e49-eee3-42b6-ba62-76c9671b9f65
# ╠═84449c10-438a-44b8-ac1e-dd617296dc8a
# ╠═e4eb0ee4-76eb-4b63-ae81-6efaf90c5338
# ╠═36dd5725-1589-4f29-ae6c-9f4ebe223e83
# ╠═7379277e-aec8-4bc6-b763-cadd0ebd6314
# ╠═a4b5f824-e256-40f8-8d59-02605f0954f7
# ╠═4c528352-37d6-4cdb-8837-b4234b526f68
# ╠═b9ebfa65-f18a-47c5-a4e6-4af7305d9261
# ╠═a886d29e-bce2-4316-860d-3063df707f36
# ╠═fabb97a8-701c-4d37-b679-a254b95b0a26
# ╠═015073d4-9e93-40df-b69b-3fd7b0026c72
# ╠═ac0da309-8776-45b6-84e3-cdd8f574ac7c
# ╠═c587da33-bd3a-4ed2-8cdc-61cbcb7b86d9
# ╠═075e3263-a63f-4172-aacb-ed47ef4b1e79
# ╠═101dcdc7-b40b-4f29-aafd-451ecfea8ac7
# ╠═01059502-8b67-41be-9a89-fd2fbb938267
# ╠═4447861d-56f6-44ae-91d6-87a393be0a43
# ╠═eb0c5ef2-28d5-431f-8d22-7db0564c24a3
# ╠═6a6ec27c-4cd1-429b-b881-f92552550118
# ╠═76566bd0-00b0-4ce4-871a-ed6478dd03f1
# ╠═3a70587b-4902-471c-8516-5fa79d2c3e4f
# ╠═2278063b-9608-4215-8d7e-fbeb948ae6dc
# ╠═0ad4b081-015c-4544-8db0-1e1f2e6f8c22
# ╠═63e2463d-fa27-44cf-a2cc-fc53caa019f6
# ╠═13b14998-2bf9-4e5d-b902-d480a5d9711c
# ╠═f1cda707-f1c6-479c-80a7-fdf13e818528
# ╠═9168fe26-f50d-42c5-9947-8b7af9701257
# ╠═dd6ebd84-d654-4b85-94cf-5f08be08e2ac
# ╠═2b9ea586-69c2-41b8-97f3-2be58186652b
# ╠═24ddccd5-8787-4c58-8c62-4dca6787762f
# ╠═83476f5b-9e21-40e0-a91f-32f07e42af9a
# ╠═422ba08d-fea2-41fe-b599-0ba01fe1f0a6
# ╠═03399151-21e9-4e5f-aa6e-a7813ca47a15
# ╠═cb514a8a-123c-4839-a17d-8d1fb41430ab
# ╠═dd6451df-c02e-4511-ab39-fcf5c1075cd2
# ╠═59d66bf5-20d9-4eba-b245-28341aec58c0
# ╠═0be7d6ab-313c-4251-a3fc-014d0ff08d8e
# ╠═a83b79a2-176a-4055-9e85-c7ac48b80426
# ╠═fa662437-2904-4eaa-b686-ad690d5775e0
# ╠═50d8b353-12d0-45f2-833a-0199f18b7b1f
# ╠═90f6ddd0-76b8-442e-a976-84a0f8233c7c
# ╠═8baed9ad-8824-41f8-ab8a-c2acf60edb46
# ╠═0a8fad3e-759e-4b8b-a3a0-466574feb9c9
# ╠═4f527a3c-5e0d-4b7d-916f-ff580eb8749e
# ╠═1cd02f96-6e4e-44d9-b44d-a15147b626ba
# ╠═fad5b501-185b-4a0a-863c-2fcebb933385
# ╠═6df45b5a-69b1-45d9-985d-72cd6cd4f55e
# ╠═f8653bd3-d2f0-43fd-9d64-381aa7a5c5aa
# ╠═dfe40c17-4f5d-4a63-b6d1-cb4d6f4d178d
# ╠═9254e969-0f77-4cdb-9198-618fee92d1db
