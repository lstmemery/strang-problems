### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 391180f4-a3a0-11eb-19d9-9b55159b8eb8
using LinearAlgebra

# ╔═╡ 862350fa-41d1-43ba-84f3-1fdbc439a5da
using Symbolics

# ╔═╡ 0851a8f3-344f-4917-b90b-3788f559e21f
using Latexify

# ╔═╡ 504449e8-d42d-4bc9-aae8-6b3864cf682f
using LaTeXStrings

# ╔═╡ 08eb66bb-62eb-4062-8131-f24f86e43a0b
using Permutations

# ╔═╡ 80669fad-ed01-40c2-9b13-a1fc0c456edb
using SpecialMatrices

# ╔═╡ cdf51b0c-96f7-400d-8032-18a660f71482
function determinant(A)
	A_lu = lu(convert(Matrix, A), check=false)
	sign(Permutation(A_lu.p)) * reduce(*, diag(A_lu.U))
end

# ╔═╡ cd6fbb72-5698-4131-b356-0e368866eac1
two_by_two_det(A) = A[1,1]*A[2,2] - A[1,2]*A[2,1]

# ╔═╡ bc146f81-dc5f-4f4d-91a5-4cf781084484
rat(a) = rationalize(a, tol=1e-4)

# ╔═╡ 167d7783-d869-4bd6-874f-17f7a17f2702
md"""
``\newcommand\m[1]{\begin{bmatrix}#1\end{bmatrix}}``
"""

# ╔═╡ 0122ca3c-a230-427f-960e-c62e514f9a3c
md"""
## Q1

``\det(A) = \frac{1}{2}``

``\det(2A) = \frac{2^4}{2} = 8``

``\det(-A) = -\frac{1^4}{2} = \frac{1}{2}``

``\det(A^2) = \det(A)\det(A) = \frac{1}{4}``

``\det(A^{-1}) = 2``
"""

# ╔═╡ 0df04316-405d-401f-af68-7ee8bb85ef9d
md"""
## Q2

``\det(A) = -1``

``\det(\frac{1}{2}A) = (\frac{-1}{2})^3 = -\frac{1}{8}``

``\det(-A) = 1^3 = 1``

``\det(A^2) = 1``

``\det(A^{-1}) = -1``
"""

# ╔═╡ 4d2b5dee-fd01-4789-bb56-6b0026cf9bd3
md"""
## Q3
"""

# ╔═╡ 1791e6c6-ae75-4172-9c83-14f882b38cde
# False
q2AI = I(2)

# ╔═╡ ebaf7468-96f9-44cc-b36a-c1e1657010e3
q2AA = [
	1 0;
	0 2
]

# ╔═╡ 7ea5de03-abb6-49d4-88f6-eedbf57d8877
determinant(q2AI + q2AA) == 1 + determinant(q2AA)

# ╔═╡ 68bf5c07-8dca-487f-826c-83dcb15378fc
md"""
b) True because ``\det(AB) = \det(A)\det(B)``. After ``\det(AB)`` is found we could simply rename ``AB`` to ``A`` and ``C`` to ``B`` get the final determinant.
"""

# ╔═╡ e6bf8675-7e02-4faf-8f40-4bfad675f1fd
# c) False
q2CA = 2 * I(2)

# ╔═╡ 4320e1bc-9d37-40cf-9b35-7801108cee0f
determinant(4*q2CA)

# ╔═╡ 29763ce5-4bc2-46af-b1ae-12c9ec21094b
4*determinant(q2CA)

# ╔═╡ c39eb249-c07c-4ad5-a88f-9c22daa7dc6c
q3DA = [
	0 0;
	0 1
]

# ╔═╡ 719c2529-cd8b-4dad-907f-2ed14ec204c2
q3DB = [
	1 3;
	5 7
]

# ╔═╡ 042b6dc6-f33d-4a40-b012-ee8104230e8e
det(q3DA*q3DB - q3DB*q3DA)

# ╔═╡ 57d51187-01ef-443a-a773-1eacde026ba9
md"""
## Q4
"""

# ╔═╡ d5e35a0e-30d1-4c7e-b58f-874181f30931
q4AJ3 = [
	0 0 1;
	0 1 0;
	1 0 0
]

# ╔═╡ b29d71b1-2ec6-456f-a4ee-2733179edb86
q4Ap13 = q4AJ3

# ╔═╡ 9a3ebbbc-ea20-47f5-85a3-2a442bf04b84
q4Ap13 * q4AJ3

# ╔═╡ 45bccc5b-8448-401b-adb0-ff149f11180a
q4BJ4 = [
	0 0 0 1
	0 0 1 0;
	0 1 0 0;
	1 0 0 0
]

# ╔═╡ f155ae76-5c5e-435f-9803-366a35e93cdd
q2Bp14 = [
	0 0 0 1;
	0 1 0 0;
	0 0 1 0;
	1 0 0 0
]

# ╔═╡ 57df4552-d3c6-49d2-afa5-91cff7ad74d4
q2Bp23 = [
	1 0 0 0;
	0 0 1 0;
	0 1 0 0;
	0 0 0 1
]

# ╔═╡ fb50935b-9d9c-4b42-94d9-e45d5a384c75
q2Bp14 * q2Bp23 * q4BJ4

# ╔═╡ 84941a90-ce98-42e6-8a42-18ed52ace1ae
md"""
## Q5
"""

# ╔═╡ cca7a963-0640-40a1-8960-1bb300d843d6
q5a = [determinant(reverse(convert(Matrix, I(i)), dims=1)) for i=2:101]

# ╔═╡ a4a8c363-62c7-4b8f-ab2f-68f1d894d17a
q5b = [i % 4 < 2 ? 1 : -1 for i=2:101]

# ╔═╡ f7738c8a-f65c-4890-96fd-93ed0906cc53
q5a == q5b

# ╔═╡ 2de08088-2619-45b7-87c8-cbf7c92c4a5f
md"""
## Q6

We can multiply the 0-row by any constant ``c`` but because ``c*0=0`` the determinant would not change. The only possible value that can be mutliplied by an arbitrary constant and remain the same is 0. Thus, the determinant much be 0.

Alternative:

``\det\left(\m{a - a & b - b \\ c & d}\right) = \det\left(\m{a & b \\ c & d}\right) - \det\left(\m{a & b \\ c & d}\right) = 0``
"""

# ╔═╡ b0fc40fb-3f87-40e0-b3c8-8c5d20057499
md"""
## Q7
"""

# ╔═╡ e49235ba-f400-4fec-8858-75671c122b4f
begin
	@variables θ

	q7Q1 = [
		cos(θ) -sin(θ);
		sin(θ) cos(θ)
	]
	L"%$(q7Q1)"
end

# ╔═╡ dfd280d9-a84f-4b4d-9c28-158802629481
L"%$(simplify.(determinant(q7Q1)))"

# ╔═╡ 2d45885e-97f4-4032-ade5-29557da26443
md"""
``\cos^2(\theta) + \sin^2(\theta) = 1``
"""

# ╔═╡ 44f24ca5-211d-4892-b818-586f2f0c3003
begin
	q7Q2 = [
		(1 - 2*cos(θ)^2) -2*cos(θ)*sin(θ);
		-2*cos(θ)*sin(θ) (1 - 2*sin(θ)^2)
	]

	L"%$(q7Q2)"
end

# ╔═╡ 22e40af6-d1c4-40da-af82-fc8a0937b37d
md"""
``(1 - (2(\cos(θ)^2)))*(1 - (2(\sin(θ)^2)) - (4(\cos(θ)^2)*(\sin(θ)^2)*((1 - (2(\cos(θ)^2)))^-1)))``

``1 -2\cos^2(\theta) - 2\sin^2(\theta) - 4\cos^2(\theta)\sin^2(\theta) + 4\cos^2(\theta)\sin^2(\theta)``

``1 -2(\cos^2(\theta) +\sin^2(\theta)) = 1 - 2 = -1``
"""

# ╔═╡ 99704b34-df57-4875-8473-e601f93a0da0
md"""
## Q8

a)

``Q^TQ=I``


``\det(Q^TQ) = \det(Q^2) = (\pm 1)^2 ``

b)

Since ``Q^TQ = I = Q^2`` every second power should reset the the determinant back to ``1``
"""

# ╔═╡ 615cc04a-0d13-46db-8ac3-aa3ef5b57649
md"""
## Q9
"""

# ╔═╡ 75efa487-6ea1-40e9-811d-0ab1d73fe73c
q9A = [
	0 0 1;
	1 0 0;
	0 1 0
]

# ╔═╡ 1f004966-c456-49fb-9aec-9a065a6e4d6e
determinant(q9A)

# ╔═╡ a2585148-5de4-46ab-bcc8-f88920b542ff
q9B = [
	0 1 1;
	1 0 1;
	1 1 0
]

# ╔═╡ a1d7bf3a-5e0c-4700-b9b4-c771e4772ec8
determinant(q9B)

# ╔═╡ 8975907a-a1ae-42c4-b84c-e995c254a6ac
q9C = [
	1 1 1;
	1 1 1;
	1 1 1
]

# ╔═╡ 0d6dad8c-401b-4e71-b2d7-ef6bd6a99eba
determinant(q9C)

# ╔═╡ cad26fb6-dea4-44d9-8c77-bd89777de379
md"""
## Q10
"""

# ╔═╡ 6743660c-a556-41a6-a4f3-e7135ee2d394
q10A = [
	101 -100;
	-100 101
]

# ╔═╡ b131f193-9e2b-4b78-ae47-d6195302612d
determinant(q10A - I)

# ╔═╡ b6d49975-23cb-41a1-bfc2-df6e3495a9ca
rat.(determinant(q10A))

# ╔═╡ f625babb-c2aa-4869-9a48-70bcd8f5c17d
md"""
## Q11

``CD = -DC``

The ``-1`` gets multiplied for each row of the matrix ``n``

``\det(CD) = (-1)^n\det(CD)`` when ``n \mod 2 = 0``
"""

# ╔═╡ 82824596-25b2-4f76-9b53-22c55dd009f0
md"""
## Q12

The problem with this caluculation is that the factor outside of the determinant ``\frac{1}{ad-bc}`` needs to be applied to each row. So:

``\det(A^{-1}) = \frac{ad-bc}{(ad-bc)^2} = \frac{1}{ad-bc}``
"""

# ╔═╡ 1d04dbe6-4349-4958-b767-593fb0306673
md"""
## Q13
"""

# ╔═╡ 5e48145d-dfcb-4936-a40f-d4db4e387039
q13A1 = [
	1 1 1;
	1 2 2;
	1 2 3
]

# ╔═╡ 97222571-4c7f-4d25-ab2d-7169ced54349
determinant(q13A1)

# ╔═╡ d84e412b-b530-49c3-b3d4-64f58ed8d0b0
q13A2 = [
	1 2 3;
	2 2 3;
	3 3 3
]

# ╔═╡ 6e70cb7d-8a08-4649-a4b2-a6697a972a70
determinant(q13A2)

# ╔═╡ d4ec001d-8520-410d-9957-69fc2a29bac8
md"""
## Q14
"""

# ╔═╡ 30f6fe29-6cd9-413b-872b-826c393aa00d
q14A1 = [
	1 2 3 0;
	2 6 6 1;
	-1 0 0 3;
	0 2 0 7
]

# ╔═╡ 46f6fca0-db76-4665-b05a-2a3859949ba4
determinant(q14A1)

# ╔═╡ 87df251d-ced5-485e-980a-0935473b12cf
q14A2 = [
	2 -1 0 0;
	-1 2 -1 0;
	0 -1 2 -1;
	0 0 -1 2
]

# ╔═╡ 6c383aee-fc30-4390-85c8-299676441e90
determinant(q14A2)

# ╔═╡ 48769477-4120-49bf-ad48-099e9a52d64b
md"""
## Q15
"""

# ╔═╡ 1a4eef3c-1d83-495d-8ac0-5c44aa5d5bdc
q15A1 = [
	101 201 301;
	102 202 302;
	103 203 303
]

# ╔═╡ 436e602d-7dc8-4ec5-9831-e53294e58fe9
determinant(q15A1)

# ╔═╡ 0f6a6ba0-269a-4e58-b99f-8e302cf7b234
begin
	@variables t
	q15A2 = [
	1 t t^2;
	t 1 t;
	t^2 t 1
]
	L"%$(q15A2)"
end

# ╔═╡ 0b229cd9-8449-4104-b585-d9709da0d67e
L"%$(simplify.(determinant(q15A2)))"

# ╔═╡ 4ec93e05-d6f5-423b-905b-a1f4db1a1a7c
md"""
``(1 - t^2) \left(1 - t^4 - \frac{(t - t^3)^2}{1 - t^2} \right)``

``(1 - t^2)(1 - t^4) - (t-t^3)^2``

``1 - t^2 - t^4 + t^6 - t^2 + 2t^4 - t^6 = 1 - 2t^2 + t^4 = (t^2-1)``
"""

# ╔═╡ 6fa470b1-8c85-426f-a490-671a55e22475
md"""
## Q16
"""

# ╔═╡ d3930802-1772-40db-99ab-e6d6ec4abd39
q16A1 = [1 2 3]'*[1 -4 5]

# ╔═╡ d6e6ce13-9dae-4a45-b132-2797ff9e42f0
determinant(q16A1)

# ╔═╡ 493b3421-41c8-468f-902b-25752fe1db2a
q16A2 = [
	0 1 3;
	-1 0 4;
	-3 -4 0
]

# ╔═╡ abefbe08-329f-4855-b47b-dcaa932d6239
determinant(q16A2)

# ╔═╡ 4335a4fa-3256-4583-8aaa-3f6da3394fe5
md"""
## Q17
"""

# ╔═╡ 8c0d1881-5de9-40ee-8adb-8e5d8edc2985
begin
	@variables a b c
	q17A = [
		0 a b;
		-a 0 c;
		-b -c 0
	]
	L"%$(q17A)"
end

# ╔═╡ 8032bc21-c8d1-4597-a095-0cd0882e65ee
L"%$(simplify.(determinant(q17A)))"

# ╔═╡ ce9d63d2-b1fd-467c-a227-759622e2491b
q17A2 = [
	0 2 2 1;
	-2 0 3 2;
	-2 -3 0 1;
	-1 -2 -1 0
]

# ╔═╡ 0f35fb3f-722d-4a10-8151-7f83a14c1c9d
determinant(q17A2)

# ╔═╡ 2ae7e7a9-c6c0-454a-ba39-4441421a89cd
md"""
## Q18
"""

# ╔═╡ c3230af5-972e-415a-82ee-b12ed5c01e99
begin
	q18A = [
		1 a a^2;
		1 b b^2;
		1 c c^2
	]
	L"%$(q18A)"
end

# ╔═╡ bd1916c9-327e-4626-b3a0-f4ad662e2b7b
L"""%$(lu(q18A).U)"""

# ╔═╡ 29b0d257-3060-4aa9-9f84-d2ac7ee8509d
md"""
``\m{1 & a & a^2 \\ 0 & b - a & b^2 - (a^2) \\ 0 & 0 & c^2 - (a^2) - ((c - a)*((b - a)^-1)*(b^2 - (a^2)}``


``(b -a)(c^2 - a^2 - \frac{c-a}{b-a}*(b^2 - a^2)``

``(b -a)c^2 - (b -a)a^2 - (c-a)(b^2 - a^2)``

``(b -a)c^2 - (b -a)a^2 - (c-a)(b-a)(b+a)``

``(b -a)(c^2 - a^2 - (c-a)(b+a))``

``(b -a)((c-a)(c+a) - (c-a)(b+a))``

``(b -a)(c -a)(c+a - b -a))``

``(b -a)(c -a)(c - b)``

"""

# ╔═╡ 6fedf279-d148-41b4-a416-69a54414c1f4
md"""
## Q19
"""

# ╔═╡ 08e456b2-c964-47d9-a692-39fdc9202937
q19U1 = [
	1 4 6;
	0 2 5;
	0 0 3
]

# ╔═╡ f9c538c3-76ba-46e5-b4b7-8becfe232ad2
determinant(q19U1)

# ╔═╡ e5380b20-32ad-4627-9de0-df0a3b593417
rat.(determinant(q19U1 \ I))

# ╔═╡ 4c78426d-d211-47d5-ac9a-8017440db59e
determinant(q19U1^2)

# ╔═╡ 965c4deb-3556-452a-9c27-ff8d27483486
begin 
	@variables d
	q19U2 = [
		a b;
		0 d
	]
	L"%$(q19U2)"
end

# ╔═╡ 2b64623b-3dd0-472f-81b2-57a4445f97a2
L"%$(determinant(q19U2))"

# ╔═╡ 8ff678bc-71a0-4c6e-8631-de50348ecd7a
L"%$(determinant(q19U2 \ I))"

# ╔═╡ 46f8ca33-20d6-4ac2-87aa-06a3ffca5114
L"%$(determinant(q19U2^2))"

# ╔═╡ b2af882f-eaaf-48a2-ab72-a8c3317f7800
md"""
## Q20
"""

# ╔═╡ fa42fd8b-3dee-42b2-b607-b029ea8ecc7e
begin
	@variables L l
	q20A = [
		(a - L*c) (b -L*d);
		(c - l*a) (d - l*b)
	]
	L"%$(q20A)"
end

# ╔═╡ 29033d22-7438-44c7-9949-68b2fea84c14
L"%$(determinant(q20A))"

# ╔═╡ f880ae4d-388b-48cd-9638-2a88667dcce8
md"""
``(a - Lc)\left((d - bl) - \frac{(b - Ld)(c - al)}{a - Lc}\right)``


``(a - Lc)(d - bl) - (b - Ld)(c - al)``

``ad - abl - cdL + bclL - bc + abl + cdL - adlL``

``(1 - lL)(ad - cd)``
"""

# ╔═╡ c75468e5-e757-4427-935b-ab2177f25074
md"""
## Q21
"""

# ╔═╡ 5dfff480-16ab-4176-95a0-4cfdd8c92595
md"""
## Q22
"""

# ╔═╡ ef7c5a46-8577-40de-9fb4-fc3d18eed794
q22A1 = [
	2 1;
	1 2
]

# ╔═╡ c27e3621-fdaf-4999-bd01-76c81a4807fb
two_by_two_det(q22A1)

# ╔═╡ 7f6c993b-440d-464d-989f-6a850257da91
q22A2 = (1//3) * [
	2 -1;
	-1 2
]

# ╔═╡ 9aa4a0c2-779d-4ca9-855b-483aa1db777f
two_by_two_det(q22A2)

# ╔═╡ 5786d83e-2422-4488-95ae-3ca35bdc506f
begin
	@variables λ
	q22A3 = [
		(2 -λ) 1;
		1 (2 - λ)
	]
	L"%$(q22A3)"
end

# ╔═╡ 912411a6-9151-454b-9bfa-1d4c54f9fdcf
L"%$(two_by_two_det(q22A3))"

# ╔═╡ 0d9a072e-b162-498e-8b95-70d48265dc04
Symbolics.value.(substitute.(two_by_two_det(q22A3), (Dict(λ => 1),)))

# ╔═╡ 8cf52913-65be-4f4e-8399-a29e236cde35
Symbolics.value.(substitute.(two_by_two_det(q22A3), (Dict(λ => 3),)))

# ╔═╡ 984aff93-a713-41e6-b7ff-7cfe667b73ea
md"""
## Q23
"""

# ╔═╡ efba08dc-d7ff-49a7-9ca5-6026f25fa17c
q23A = [
	4 1;
	2 3
]

# ╔═╡ 94993f4e-c956-4164-8688-de3b70bf996a
determinant(q23A^2)

# ╔═╡ bef79362-ff83-4c4d-b7e5-b9059dc6ca1e
determinant(q23A \ I)

# ╔═╡ 01363056-8824-4265-a838-ab9849f6b4c1
begin
	q23det = determinant(q23A - λ * I)
	L"%$(q23det)"
end

# ╔═╡ 70e6dd58-29e3-4659-bbb7-cf3bc033336e
md"""
``(4 - λ)(3 - λ - \frac{2}{4 - λ}) = 0``

``(4 - λ)(3 - λ) - 2 = 0``

``λ \in \{2, 5\}``
"""

# ╔═╡ 1bc67e05-7afc-44e0-9b48-942c1bb2a97b
md"""
## Q24
"""

# ╔═╡ 077af35c-cf34-4a1a-a94d-eb919c454af3
q24A = [
	3 3 4;
	6 8 7;
	-3 5 -9
]

# ╔═╡ 1d162322-1af5-4647-bc26-267bb5eef3d1
q24LU = lu(q24A)

# ╔═╡ 3f5edaee-85f2-4415-834e-c577a6b8961d
rat.(determinant(q24A))

# ╔═╡ dbe2a611-9e37-4fbf-acf3-cd2d0a02600a
rat.(determinant(q24LU.U))

# ╔═╡ 285f2270-8aee-4562-bb47-98a91d45a186
rat.(determinant(q24LU.L))

# ╔═╡ 79975186-f92a-4290-bddb-76a734dd971d
rat.(determinant(inv(q24LU.U) * inv(q24LU.L)))

# ╔═╡ c1119156-90d1-45bf-a449-0e88620382ed
rat.(determinant(inv(q24LU.U) * inv(q24LU.L) * q24A))

# ╔═╡ 4d6010a8-d5dc-4c8a-9633-19192e896b9f
md"""
## Q25

Determinants are found through elimination, which does not change the determinant. By definition, eliminating one row that is a constant multiple of another is going to produce a ``0`` row. Any matrix with a ``0`` row will have a determinant of ``0``, as well.

"""

# ╔═╡ cc965cd3-db89-4823-b229-25d9de6f96bc
md"""
## Q26

The difference between each row is the ``1`` vector. If ``n \gt 2`` then that means that three are multiple ways to get to the ``1`` vector and the matrix is singular. Singular matrices have a ``\det(A) = 0``.
"""

# ╔═╡ 73de7b80-8f00-49c2-a54e-4cdde2cbaaa1
md"""
## Q27
"""

# ╔═╡ eb4cb9ab-dce1-4d73-b7fd-7366459499f0
begin
	q27A = [
		0 a 0;
		0 0 b;
		c 0 0
	]
	L"%$(q27A)"
end

# ╔═╡ 1907d406-3e51-42a9-891a-b3b9a0b4e972
L"%$(determinant(q27A))"

# ╔═╡ 691b5bc3-6bce-4e0b-8d70-1a0fb6afc98c
begin
	q27B = [
		0 a 0 0;
		0 0 b 0;
		0 0 0 c;
		d 0 0 0
	]
	L"%$(q27B)"
end

# ╔═╡ dce2cd24-60be-4dc0-a3d3-2705955e4cab
L"%$(determinant(q27B))"

# ╔═╡ d1c2a420-17a8-44cf-bff9-29193d0f46dd
begin 
	q27C = [
		a a a;
		a b b;
		a b c
	]
	L"%$(q27C)"
end

# ╔═╡ 28cac7fa-0afc-4e7e-ad98-f09391a7611c
L"%$(determinant(q27C))"

# ╔═╡ 1fc75285-4cda-4781-a1fc-89dadd059018
md"""
## Q28
"""

# ╔═╡ 4c40e37e-42e2-44e5-a0b4-2ad7d65ae4bc
# b) False. doesn't account for row exchanges 

q28bA = [
	0 1;
	1 0
]

# ╔═╡ 330b24aa-84f1-40f7-9c29-af51c0e4e65e
determinant(q28bA)

# ╔═╡ dc60727f-256d-40c2-b3a8-dcd43100ef16
q28cA  = [
	1 0;
	1 0
]

# ╔═╡ 7116c596-4fbd-4e17-a9f1-ff1d5ecc0ecb
determinant(q28bA - q28cA)

# ╔═╡ 502f5d4c-074b-44dc-ad82-5955a6c18f5a
determinant(q28bA) - determinant(q28cA)

# ╔═╡ 8fd517e9-cf62-4a50-9b1e-baca61d2e76f
md"""
d) True!

``\det(AB) = \det(A)\det(B) = \det(B)\det(A) = \det(BA)``
"""

# ╔═╡ b6a557c0-501a-4f5e-938d-aad3d4e99881
md"""
## Q29

``A`` is not guaranteed to be square and non-square matrices do not have determinants.
"""

# ╔═╡ 9367d99d-bb0b-4c32-af48-34a4c9f10c25
md"""
## Q30

``\frac{\partial{f}}{\partial{a}} = \frac{a}{ad -bc}``


``\m{\frac{a}{ad-bc} &-\frac{c}{ad-bc} \\ -\frac{b}{ad-bc} & \frac{d}{ad-bc}} = \frac{1}{ad-bc}\m{a &-c \\ -b & d}``

``\det\left(\frac{1}{ad-bc}\m{a &-c \\ -b & d}\right) = \frac{ad-bc}{(ab-bc)^2} = \frac{1}{ad-bc}``
"""

# ╔═╡ dfe4de89-04d5-4fe9-8cea-ade42e75a83c
md"""
## Q31
"""

# ╔═╡ 3c9e0ee9-f34b-4763-8185-c83ad9360399
q31dets = map(x->det(big.(Hilbert(x))), 1:10)

# ╔═╡ fd0aa401-d42a-440f-ac6b-b6d7910219c4
diag(lu(Hilbert(5)).U)

# ╔═╡ d9e8b617-786b-4f16-9363-8cb82bc11a5b
md"""
## Q32
"""

# ╔═╡ 2c562b50-7efb-460c-9bb1-6fcf348214d1
rand(Float64,(50,50))

# ╔═╡ cebe6a7a-a946-40b3-a506-cc998e88be64
function det_average(gen::Function, size::Int, iterations::Int)
	reduce(+,map(x -> det(gen(Float64, (size,size))), 1:iterations))/size
end

# ╔═╡ a83c3243-7ad6-46bc-8704-d8cf027b5a25
det_average(rand, 50, 100)

# ╔═╡ 647e8408-efab-4d2c-a454-19ea92efdad3
det_average(rand, 100, 100)

# ╔═╡ fc4df60d-530e-4e9d-b7cb-3d9c2b41d0b4
det_average(rand, 200, 100)

# ╔═╡ 745b87f3-ade8-47ba-86c8-af2f6a2480da
det_average(rand, 400, 100)

# ╔═╡ 7c7b79fe-6636-40a5-91f7-c559b80dd279
det_average(randn, 50, 100)

# ╔═╡ 1b007f09-937d-49e4-9c1a-897800b30be1
det_average(randn, 100, 100)

# ╔═╡ 1a00706a-ae3d-40ac-8569-6ac64815002b
det_average(randn, 200, 100)

# ╔═╡ b5868b57-9cee-4857-bb38-0c3ae1dcb93f
det_average(randn, 400, 100)

# ╔═╡ 56c241d5-e6a4-46bf-8501-15d1f9fc47a1
md"""
## Q33
"""

# ╔═╡ a1390682-acfc-4267-bd91-c183b79460f7
reshape(map(x-> x=='0' ? -1 : 1, collect(bitstring(UInt64(123235))[29:end])), (6,6))

# ╔═╡ cb84f0b8-47b3-49b3-9764-4755c418c628
q33A = [
  1   1   1   1   1  1;
 -1   1   1   1   1  1;
 -1  -1   1   -1   1  1;
 -1  -1  -1   1   1  1;
 -1  -1  -1  -1   1  1;
 1  1  1  -1  -1  1
]

# ╔═╡ e4d1fa39-fd6f-4d73-b19b-121f94b477f5
[det(Matrix(Permutation(6, i))'*q33A) for i in 1:factorial(6)]

# ╔═╡ 426c5cd5-3a60-4518-9550-0cb3ccdb277a
reshape(map(x-> x=='0' ? -1 : 1, collect(bitstring(UInt64(123235))[64-(6^2-1):end])), (6,6))

# ╔═╡ 92ffa358-1d0f-4980-955e-f9d33ad16d14
q33size = 2

# ╔═╡ c12b73cb-533f-4804-86f4-ee25cf0179c0
function maximal_determinant(q33size)
	max_determinant = 0
	max_determinant_matrix = zeros(q33size, q33size)
	show(max_determinant)

	for q33_value = 1:2^q33size^2
		q33_matrix = reshape(map(x-> x=='0' ? -1 : 1, collect(bitstring(UInt64(q33_value))[64-(q33size^2-1):end])), (q33size,q33size))
		if det(q33_matrix) > max_determinant
			max_determinant = det(q33_matrix)
			max_determinant_matrix = q33_matrix
		end
	end
	max_determinant, max_determinant_matrix
end

# ╔═╡ 2f80b7fe-6e1b-45bc-add9-43923b72a17e
md"""
This is a rough problem. Note that negations and permutations of rows and columns do not change ``|\det(A)|``. This means any `\{1, -1\}`- matrix can be converted into normalized form, where the first row and column of the matrix are all 1s. There is an additional great simplication found here: https://en.wikipedia.org/wiki/Hadamard's_maximal_determinant_problem

That says that any normalized ``{1,-1}`` matrix corresponds to a ``{0 ,1}`` matrix of ``(n-1,n-1)``. You better believe I'm going to exploit that.
"""

# ╔═╡ 4ee63a89-924f-431e-a1e1-8a7eb1281dc2
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

# ╔═╡ cdf7e11d-af82-4f7c-bb17-da08463285bd


# ╔═╡ e41a0342-9d57-4d62-b97e-2409ec07fea6
q33_maxdet, q33_maxdet_matrix = smart_maximal_determinant(6)

# ╔═╡ de351f83-010f-45bd-b154-1f83c9708782
det(q33_maxdet_matrix)

# ╔═╡ cb67a78e-8b85-4068-b534-7a96e9278c85
md"""
## Q34
"""

# ╔═╡ d511403a-1a04-4e6c-be1d-22be09c1a18f
q34A = [
	1 0 0;
	0 2 0;
	0 0 3
]

# ╔═╡ 2c9fd1b1-fc5e-4941-a146-12ebd31d3a70
determinant(q34A)

# ╔═╡ e5c2e27b-2928-47d7-bf39-f139037dc020
q34B = [
	1 2 3;
	1 2 0;
	1 0 0
]

# ╔═╡ 7f59f2aa-579f-469a-820d-0b204efa9651
determinant(q34B)

# ╔═╡ Cell order:
# ╠═391180f4-a3a0-11eb-19d9-9b55159b8eb8
# ╠═862350fa-41d1-43ba-84f3-1fdbc439a5da
# ╠═0851a8f3-344f-4917-b90b-3788f559e21f
# ╠═504449e8-d42d-4bc9-aae8-6b3864cf682f
# ╠═08eb66bb-62eb-4062-8131-f24f86e43a0b
# ╠═80669fad-ed01-40c2-9b13-a1fc0c456edb
# ╠═cdf51b0c-96f7-400d-8032-18a660f71482
# ╠═cd6fbb72-5698-4131-b356-0e368866eac1
# ╠═bc146f81-dc5f-4f4d-91a5-4cf781084484
# ╠═167d7783-d869-4bd6-874f-17f7a17f2702
# ╠═0122ca3c-a230-427f-960e-c62e514f9a3c
# ╠═0df04316-405d-401f-af68-7ee8bb85ef9d
# ╠═4d2b5dee-fd01-4789-bb56-6b0026cf9bd3
# ╠═1791e6c6-ae75-4172-9c83-14f882b38cde
# ╠═ebaf7468-96f9-44cc-b36a-c1e1657010e3
# ╠═7ea5de03-abb6-49d4-88f6-eedbf57d8877
# ╠═68bf5c07-8dca-487f-826c-83dcb15378fc
# ╠═e6bf8675-7e02-4faf-8f40-4bfad675f1fd
# ╠═4320e1bc-9d37-40cf-9b35-7801108cee0f
# ╠═29763ce5-4bc2-46af-b1ae-12c9ec21094b
# ╠═c39eb249-c07c-4ad5-a88f-9c22daa7dc6c
# ╠═719c2529-cd8b-4dad-907f-2ed14ec204c2
# ╠═042b6dc6-f33d-4a40-b012-ee8104230e8e
# ╠═57d51187-01ef-443a-a773-1eacde026ba9
# ╠═d5e35a0e-30d1-4c7e-b58f-874181f30931
# ╠═b29d71b1-2ec6-456f-a4ee-2733179edb86
# ╠═9a3ebbbc-ea20-47f5-85a3-2a442bf04b84
# ╠═45bccc5b-8448-401b-adb0-ff149f11180a
# ╠═f155ae76-5c5e-435f-9803-366a35e93cdd
# ╠═57df4552-d3c6-49d2-afa5-91cff7ad74d4
# ╠═fb50935b-9d9c-4b42-94d9-e45d5a384c75
# ╠═84941a90-ce98-42e6-8a42-18ed52ace1ae
# ╠═cca7a963-0640-40a1-8960-1bb300d843d6
# ╠═a4a8c363-62c7-4b8f-ab2f-68f1d894d17a
# ╠═f7738c8a-f65c-4890-96fd-93ed0906cc53
# ╠═2de08088-2619-45b7-87c8-cbf7c92c4a5f
# ╠═b0fc40fb-3f87-40e0-b3c8-8c5d20057499
# ╠═e49235ba-f400-4fec-8858-75671c122b4f
# ╠═dfd280d9-a84f-4b4d-9c28-158802629481
# ╠═2d45885e-97f4-4032-ade5-29557da26443
# ╠═44f24ca5-211d-4892-b818-586f2f0c3003
# ╠═22e40af6-d1c4-40da-af82-fc8a0937b37d
# ╠═99704b34-df57-4875-8473-e601f93a0da0
# ╠═615cc04a-0d13-46db-8ac3-aa3ef5b57649
# ╠═75efa487-6ea1-40e9-811d-0ab1d73fe73c
# ╠═1f004966-c456-49fb-9aec-9a065a6e4d6e
# ╠═a2585148-5de4-46ab-bcc8-f88920b542ff
# ╠═a1d7bf3a-5e0c-4700-b9b4-c771e4772ec8
# ╠═8975907a-a1ae-42c4-b84c-e995c254a6ac
# ╠═0d6dad8c-401b-4e71-b2d7-ef6bd6a99eba
# ╠═cad26fb6-dea4-44d9-8c77-bd89777de379
# ╠═6743660c-a556-41a6-a4f3-e7135ee2d394
# ╠═b131f193-9e2b-4b78-ae47-d6195302612d
# ╠═b6d49975-23cb-41a1-bfc2-df6e3495a9ca
# ╠═f625babb-c2aa-4869-9a48-70bcd8f5c17d
# ╠═82824596-25b2-4f76-9b53-22c55dd009f0
# ╠═1d04dbe6-4349-4958-b767-593fb0306673
# ╠═5e48145d-dfcb-4936-a40f-d4db4e387039
# ╠═97222571-4c7f-4d25-ab2d-7169ced54349
# ╠═d84e412b-b530-49c3-b3d4-64f58ed8d0b0
# ╠═6e70cb7d-8a08-4649-a4b2-a6697a972a70
# ╠═d4ec001d-8520-410d-9957-69fc2a29bac8
# ╠═30f6fe29-6cd9-413b-872b-826c393aa00d
# ╠═46f6fca0-db76-4665-b05a-2a3859949ba4
# ╠═87df251d-ced5-485e-980a-0935473b12cf
# ╠═6c383aee-fc30-4390-85c8-299676441e90
# ╠═48769477-4120-49bf-ad48-099e9a52d64b
# ╠═1a4eef3c-1d83-495d-8ac0-5c44aa5d5bdc
# ╠═436e602d-7dc8-4ec5-9831-e53294e58fe9
# ╠═0f6a6ba0-269a-4e58-b99f-8e302cf7b234
# ╠═0b229cd9-8449-4104-b585-d9709da0d67e
# ╠═4ec93e05-d6f5-423b-905b-a1f4db1a1a7c
# ╠═6fa470b1-8c85-426f-a490-671a55e22475
# ╠═d3930802-1772-40db-99ab-e6d6ec4abd39
# ╠═d6e6ce13-9dae-4a45-b132-2797ff9e42f0
# ╠═493b3421-41c8-468f-902b-25752fe1db2a
# ╠═abefbe08-329f-4855-b47b-dcaa932d6239
# ╠═4335a4fa-3256-4583-8aaa-3f6da3394fe5
# ╠═8c0d1881-5de9-40ee-8adb-8e5d8edc2985
# ╠═8032bc21-c8d1-4597-a095-0cd0882e65ee
# ╠═ce9d63d2-b1fd-467c-a227-759622e2491b
# ╠═0f35fb3f-722d-4a10-8151-7f83a14c1c9d
# ╠═2ae7e7a9-c6c0-454a-ba39-4441421a89cd
# ╠═c3230af5-972e-415a-82ee-b12ed5c01e99
# ╠═bd1916c9-327e-4626-b3a0-f4ad662e2b7b
# ╠═29b0d257-3060-4aa9-9f84-d2ac7ee8509d
# ╠═6fedf279-d148-41b4-a416-69a54414c1f4
# ╠═08e456b2-c964-47d9-a692-39fdc9202937
# ╠═f9c538c3-76ba-46e5-b4b7-8becfe232ad2
# ╠═e5380b20-32ad-4627-9de0-df0a3b593417
# ╠═4c78426d-d211-47d5-ac9a-8017440db59e
# ╠═965c4deb-3556-452a-9c27-ff8d27483486
# ╠═2b64623b-3dd0-472f-81b2-57a4445f97a2
# ╠═8ff678bc-71a0-4c6e-8631-de50348ecd7a
# ╠═46f8ca33-20d6-4ac2-87aa-06a3ffca5114
# ╠═b2af882f-eaaf-48a2-ab72-a8c3317f7800
# ╠═fa42fd8b-3dee-42b2-b607-b029ea8ecc7e
# ╠═29033d22-7438-44c7-9949-68b2fea84c14
# ╠═f880ae4d-388b-48cd-9638-2a88667dcce8
# ╠═c75468e5-e757-4427-935b-ab2177f25074
# ╠═5dfff480-16ab-4176-95a0-4cfdd8c92595
# ╠═ef7c5a46-8577-40de-9fb4-fc3d18eed794
# ╠═c27e3621-fdaf-4999-bd01-76c81a4807fb
# ╠═7f6c993b-440d-464d-989f-6a850257da91
# ╠═9aa4a0c2-779d-4ca9-855b-483aa1db777f
# ╠═5786d83e-2422-4488-95ae-3ca35bdc506f
# ╠═912411a6-9151-454b-9bfa-1d4c54f9fdcf
# ╠═0d9a072e-b162-498e-8b95-70d48265dc04
# ╠═8cf52913-65be-4f4e-8399-a29e236cde35
# ╠═984aff93-a713-41e6-b7ff-7cfe667b73ea
# ╠═efba08dc-d7ff-49a7-9ca5-6026f25fa17c
# ╠═94993f4e-c956-4164-8688-de3b70bf996a
# ╠═bef79362-ff83-4c4d-b7e5-b9059dc6ca1e
# ╠═01363056-8824-4265-a838-ab9849f6b4c1
# ╠═70e6dd58-29e3-4659-bbb7-cf3bc033336e
# ╠═1bc67e05-7afc-44e0-9b48-942c1bb2a97b
# ╠═077af35c-cf34-4a1a-a94d-eb919c454af3
# ╠═1d162322-1af5-4647-bc26-267bb5eef3d1
# ╠═3f5edaee-85f2-4415-834e-c577a6b8961d
# ╠═dbe2a611-9e37-4fbf-acf3-cd2d0a02600a
# ╠═285f2270-8aee-4562-bb47-98a91d45a186
# ╠═79975186-f92a-4290-bddb-76a734dd971d
# ╠═c1119156-90d1-45bf-a449-0e88620382ed
# ╠═4d6010a8-d5dc-4c8a-9633-19192e896b9f
# ╠═cc965cd3-db89-4823-b229-25d9de6f96bc
# ╠═73de7b80-8f00-49c2-a54e-4cdde2cbaaa1
# ╠═eb4cb9ab-dce1-4d73-b7fd-7366459499f0
# ╠═1907d406-3e51-42a9-891a-b3b9a0b4e972
# ╠═691b5bc3-6bce-4e0b-8d70-1a0fb6afc98c
# ╠═dce2cd24-60be-4dc0-a3d3-2705955e4cab
# ╠═d1c2a420-17a8-44cf-bff9-29193d0f46dd
# ╠═28cac7fa-0afc-4e7e-ad98-f09391a7611c
# ╠═1fc75285-4cda-4781-a1fc-89dadd059018
# ╠═4c40e37e-42e2-44e5-a0b4-2ad7d65ae4bc
# ╠═330b24aa-84f1-40f7-9c29-af51c0e4e65e
# ╠═dc60727f-256d-40c2-b3a8-dcd43100ef16
# ╠═7116c596-4fbd-4e17-a9f1-ff1d5ecc0ecb
# ╠═502f5d4c-074b-44dc-ad82-5955a6c18f5a
# ╠═8fd517e9-cf62-4a50-9b1e-baca61d2e76f
# ╠═b6a557c0-501a-4f5e-938d-aad3d4e99881
# ╠═9367d99d-bb0b-4c32-af48-34a4c9f10c25
# ╠═dfe4de89-04d5-4fe9-8cea-ade42e75a83c
# ╠═3c9e0ee9-f34b-4763-8185-c83ad9360399
# ╠═fd0aa401-d42a-440f-ac6b-b6d7910219c4
# ╠═d9e8b617-786b-4f16-9363-8cb82bc11a5b
# ╠═2c562b50-7efb-460c-9bb1-6fcf348214d1
# ╠═cebe6a7a-a946-40b3-a506-cc998e88be64
# ╠═a83c3243-7ad6-46bc-8704-d8cf027b5a25
# ╠═647e8408-efab-4d2c-a454-19ea92efdad3
# ╠═fc4df60d-530e-4e9d-b7cb-3d9c2b41d0b4
# ╠═745b87f3-ade8-47ba-86c8-af2f6a2480da
# ╠═7c7b79fe-6636-40a5-91f7-c559b80dd279
# ╠═1b007f09-937d-49e4-9c1a-897800b30be1
# ╠═1a00706a-ae3d-40ac-8569-6ac64815002b
# ╠═b5868b57-9cee-4857-bb38-0c3ae1dcb93f
# ╠═56c241d5-e6a4-46bf-8501-15d1f9fc47a1
# ╠═a1390682-acfc-4267-bd91-c183b79460f7
# ╠═cb84f0b8-47b3-49b3-9764-4755c418c628
# ╠═e4d1fa39-fd6f-4d73-b19b-121f94b477f5
# ╠═426c5cd5-3a60-4518-9550-0cb3ccdb277a
# ╠═92ffa358-1d0f-4980-955e-f9d33ad16d14
# ╠═c12b73cb-533f-4804-86f4-ee25cf0179c0
# ╠═2f80b7fe-6e1b-45bc-add9-43923b72a17e
# ╠═4ee63a89-924f-431e-a1e1-8a7eb1281dc2
# ╠═cdf7e11d-af82-4f7c-bb17-da08463285bd
# ╠═e41a0342-9d57-4d62-b97e-2409ec07fea6
# ╠═de351f83-010f-45bd-b154-1f83c9708782
# ╠═cb67a78e-8b85-4068-b534-7a96e9278c85
# ╠═d511403a-1a04-4e6c-be1d-22be09c1a18f
# ╠═2c9fd1b1-fc5e-4941-a146-12ebd31d3a70
# ╠═e5c2e27b-2928-47d7-bf39-f139037dc020
# ╠═7f59f2aa-579f-469a-820d-0b204efa9651
