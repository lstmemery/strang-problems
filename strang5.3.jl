### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 26dfe76b-fbe3-4116-8129-a5538615df66
using LinearAlgebra

# ╔═╡ 20eb19f5-d2b1-425b-bd2e-4d142d28cd60
using LaTeXStrings

# ╔═╡ 923a3522-b6f8-434d-8bb9-85d8357df51d
using Symbolics

# ╔═╡ 1644d02b-a035-4cc3-870e-5913b6bb3489
using SymbolicUtils

# ╔═╡ 763574f2-9d7c-441d-9fae-3d1253e39db7
using SymbolicUtils.Rewriters

# ╔═╡ 15495da5-b616-4540-8d60-6251ebcd8575
using Plots

# ╔═╡ 39d59b60-bcc9-4c30-a08d-753913e68d91
using Hadamard

# ╔═╡ 82ec58dd-0373-4438-86e4-4144bb0b0b34
show_num(i) = L"%$(latexify(i))"

# ╔═╡ dc93074a-ebfe-4f4e-88f4-485f5c173d43
rat(a) = rationalize(a, tol=1e-4)

# ╔═╡ 43cc9465-3f44-4aad-a45d-5f1da3e67813
md"""
``\newcommand\m[1]{\begin{bmatrix}#1\end{bmatrix}}``
"""

# ╔═╡ 58195e48-1e8a-4d16-bdc3-736ce4920a66
md"""
``\newcommand\v[1]{\begin{vmatrix}#1\end{vmatrix}}``
"""

# ╔═╡ fc2d2f95-edc8-41e6-8a10-b1173ddb51b0
cofactor_matrix(A, i, j) = A[1:end .!= i, 1:end .!= j]

# ╔═╡ b207a078-2bfc-416b-afa2-6ae646647beb
function cofactor_determinant(A)
	if size(A)[1] == 1
		return A[1, 1]
	else
		return sum([A[1,j]*(((-1)^(1 + j)) * cofactor_determinant(cofactor_matrix(A,1,j))) for j=1:size(A)[2]])
	end
end

# ╔═╡ 94599138-d357-43f1-876d-1073c0a49714
first_row_cofactors(A) = [A[1,j]*(((-1)^(1 + j)) * cofactor_determinant(cofactor_matrix(A,1,j))) for j=1:size(A)[2]]

# ╔═╡ 5742bab2-f9fb-415b-a493-6e28b46d257b
all_cofactors(A) = [(((-1)^(i + j)) * cofactor_determinant(cofactor_matrix(A,i,j))) for i=1:size(A)[1], j=1:size(A)[2]]

# ╔═╡ 6a1227cb-3248-44ef-8ae1-d1842fb3d22f
cramer(A::Matrix) = all_cofactors(A)' / det(A)

# ╔═╡ ab77210b-e692-48e6-a194-7b5562fa2e13
function bmatrix(A, b, i)
	B = copy(A)
	B[:,i]= b
	B
end

# ╔═╡ 4f1d0cb1-57a0-4dff-996c-372558de9aef
function cramer(A::Matrix, b)
	[cofactor_determinant(bmatrix(A,b,i))/cofactor_determinant(A) for i=1:size(b)[1]]
end

# ╔═╡ fc5aef5a-1e0b-4316-86db-af2ab78110dc
begin
	@variables x a b c d e f g h i j k l m n o p q r s t u v w y z D θ ϕ ρ J u1 u2 u3 v1 v2 v3 w1 w2 w3
end;

# ╔═╡ 98422293-89b8-4f80-830b-9c360f52c42d
function cross_product(u, v)
	A = vcat([i j k], u,v)
	cofactor_determinant(A)
end

# ╔═╡ c2b022d8-2ecf-43cf-844d-4b9e4bc198c4
function cross_product_non_symbolic(u, v)
	A = vcat([1 1 1], u,v)
	first_row_cofactors(A)
end

# ╔═╡ c0b1135e-90b3-4308-b668-4ddd3c5b97fe
cross_product_non_symbolic([3 2 0], [1 4 0])

# ╔═╡ c6ff16a6-96b8-4f45-af2d-05aca7d97d84
cross_product_non_symbolic([1 1 1], [1 1 2])

# ╔═╡ 3272d2ca-50c0-4c2f-a363-935af2bf00b3
function cross_product_length(u,v)
	cross = cross_product(u,v)
	sqrt(sum(values(Symbolics.value(cross).dict) .^ 2))
end

# ╔═╡ 52a08df7-5625-4833-9ab7-793534697114
function triangular_area(u, v, w)
	tri_matrix = hvcat((3, 3, 3), u..., 1, v..., 1, w..., 1)
	abs(cofactor_determinant(tri_matrix) / 2)
end

# ╔═╡ 70bd870f-cc32-4840-8ea1-f54bb02855e4
function two_by_two_inv(A)
	a,b,c,d = A[1,1], A[1, 2], A[2, 1], A[2, 2]
	(1/(a*d - b*c)) .* [d -b; -c a]
end

# ╔═╡ 56fc6a16-5502-4116-a507-1559b6672786
vector_length(v) = sqrt(reduce(+, v .^2))

# ╔═╡ 1b3b179c-2aa2-4fa4-9eaa-0ccfcd64ed4b
triple_product(u, v, w) = (cross_product_non_symbolic(u, v)' * w')[1]

# ╔═╡ e2993de4-bd70-11eb-0ad0-75d1f93d8a80
md"""
## Q1
"""

# ╔═╡ 0f017eb0-a10e-420a-9d68-4a77348b1c26
q1aA =[
	2 5;
	1 4
]

# ╔═╡ 802f6c24-5657-4e27-bfab-d340ee6111ff
q1ab = [1 2]'

# ╔═╡ c41c25c3-3b62-43e1-aada-b90c3de9efed
cramer(q1aA, q1ab)

# ╔═╡ 9e6f0924-ce5e-4827-a4f9-d8f47d2f2fea
q1bA = [
	2 1 0;
	1 2 1;
	0 1 2
]

# ╔═╡ c064ad99-9ce2-4981-8607-04c1c7a5d85a
q1bb = [1 0 0]'

# ╔═╡ ca6821e1-bff3-4aa6-b9c7-922f258ae55c
rat.(cramer(q1bA, q1bb))

# ╔═╡ 9660772f-bbe6-4e63-90d5-212862840c2b
md"""
## Q2
"""

# ╔═╡ f3bc7f53-00bb-4ede-b353-189fd41c60fd
q2aA = [
	a b;
	c d
]

# ╔═╡ b54a70e9-a45f-4156-b51b-b8fb4e030c32
q2ab = [1 0]'

# ╔═╡ d230dc4e-6ef2-4919-95ed-9f2603bcfd9f
det(bmatrix(q2aA,q2ab,2))/det(q2aA)

# ╔═╡ 2e725c38-318c-479b-901e-98733856df60
q2bA = [
	a b c;
	d e f;
	g h i
]

# ╔═╡ ecdcac60-2410-4038-af98-a51d67f3e266
q2bb = [1 0 0]'

# ╔═╡ 1cf8bc16-a088-4f82-b414-5c381376e744
substitute(cofactor_determinant(bmatrix(q2bA,q2bb,2))/det(q2bA), Dict([det(q2bA) => D]))

# ╔═╡ 56b0a487-a01a-42dd-9382-db38c5d4e001
det(bmatrix(q2bA,q2bb,2))

# ╔═╡ e2344c48-5145-45d3-a26b-7128e20f9316
md"""
## Q3
"""

# ╔═╡ 4ff64552-e76d-4d8f-89d1-d8b3d0f98aae
q3aA = [
	2 3;
	4 6
]

# ╔═╡ da6ed2e3-0a1a-4126-81a8-2f077823ce87
cofactor_determinant(q3aA)

# ╔═╡ 1b2b7814-6550-454c-b323-8b776991940a
q3ab = [1 1]'

# ╔═╡ 3e2c4dc5-87b4-479d-9a79-f7e986cdf412
cofactor_determinant(bmatrix(q3aA,q3ab,1)) # x_{ij} / 0

# ╔═╡ ffda5ac4-0964-45f0-b0ed-45b99c5cc883
q3bb = [1 2]'

# ╔═╡ 1ead8e56-5791-4951-ad89-d26ee94e3bec
cofactor_determinant(bmatrix(q3aA,q3bb,1)) 

# ╔═╡ c4f284af-c4ee-4e3b-84ac-e2d01d2a2f58
md"""
## Q4

a)
``x_1 = \frac{\det(B_1)}{\det(A)} = \frac{|b \\ a_2 \\ a_3|}{\det(A)}``

b) Substitute ``b \rightarrow x_1a_1 + x_2a_2 + x_3a_3``

Because this operation is linear this gives us:

``\v{x_1a_1 + x_2a_2 + x_3a_3 & a_2 & a_3} = x_1\v{a_1 & a_2 & a_3} + x_2\v{a_2 & a_2 & a_3} + x_3\v{a_3 & a_2 & a_3} = x_1\v{a_1 & a_2 & a_3} + 0 + 0 ``

"""

# ╔═╡ 7f5612c7-1e5a-4206-8f0b-927b78ba6625
md"""
## Q5
"""

# ╔═╡ 9f3a3a10-d447-49ae-9ed2-22a0d96235a5
q5A = [
	a b c;
	d e f;
	g h i
]

# ╔═╡ 96387e53-dc09-403a-b900-7604ff6ecbe3
q5b = [a d g]'

# ╔═╡ 18fa8409-ac1c-4d56-b6d9-26f0ad30214b
q5_cramer = cramer(q5A, q5b)

# ╔═╡ f0044b20-a546-4c88-8f0b-6a16c98ac33a
md"""
Notice that the numerator of ``x_2`` and ``x_3`` is zero! It's the same column twice. It also makes sense that ``x_1 = 1`` because ``B_1 = A`` so ``x_1 = \frac{|B_1|}{|A|} = 1``
"""

# ╔═╡ 48e7e5d3-305b-4591-9cab-d8b70ee9b595
simplify.(q5_cramer)

# ╔═╡ 483c1ff0-85f2-469e-ba15-80910ba80fd4
md"""
## Q6
"""

# ╔═╡ 8861908e-0442-426b-940f-3b4ba9996a18
q6aA = [
	1 2 0;
	0 3 0;
	0 7 1
]

# ╔═╡ bcd49c46-388a-4d64-846e-1f81d38c2f6a
rat.(cramer(q6aA))

# ╔═╡ 6d661775-d16d-4319-bfc7-c905b5fcc242
rat.(q6aA\I)

# ╔═╡ d72684f0-7dbf-4875-b974-f41b6f530b2f
q6bA = [
	2 -1 0;
	-1 2 -1;
	0 -1 2
]

# ╔═╡ 0e64596d-8408-4d36-b7c8-41bcca851bb9
rat.(cramer(q6bA))

# ╔═╡ c6727867-89a4-4cb6-9858-eb4acea7b338
rat.(q6bA\I)

# ╔═╡ ead0878d-fcea-4666-ba11-b665e7a49880
md"""
## Q7

If all cofactors are 0 then ``C^T = 0``. Having all non-zero cofactors does not guarantee that the matrix has an inverse. There could still be no inverse if ``det(A) = 0``.
"""

# ╔═╡ 29fc5286-003d-41f1-9b68-404ff0537c40
md"""
## Q8
"""

# ╔═╡ fbddb9b0-ac15-4fee-bb18-c15d2951a083
q8A = [
	1 1 4;
	1 2 2;
	1 2 5
]

# ╔═╡ 7bbf2f1d-c963-4f9d-ae88-ea0fa80f10e8
q8C = all_cofactors(q8A)

# ╔═╡ dc319874-3405-4a50-a7b6-a93dd9b3ca0a
q8A * q8C'

# ╔═╡ 962a0d8a-bff8-49ac-bf56-5aec93f154ff
q8Ab = [
	1 1 100;
	1 2 2;
	1 2 5
]

# ╔═╡ 5835b954-24d0-4d2c-84af-e020796c5e5f
q8Cb = all_cofactors(q8Ab)

# ╔═╡ 52efbead-2202-4fd3-b3da-1fa9580fbe16
q8Ab * q8Cb'

# ╔═╡ 66041ea1-6876-47f8-b3db-42a722b4531b
det(q8A)

# ╔═╡ 30939bc8-6eca-48bd-8e11-bb68de9428fc
md"""
## Q9

When ``\det(A) = 1 \rightarrow A^{-1} = C^T \rightarrow A = C^{-T}``
"""

# ╔═╡ a2a52c1e-ec88-488e-b1c4-185f15a0bea5
md"""
## Q10


``AC^T = (\det(A))I``

``\det(AC^T) = \det(A)\det(C^T) = \det(A)^n``

``\det(C^T) = \det(C) = \det(A)^{n-1}``
"""

# ╔═╡ 19f93fb4-6b2a-4c19-b755-6ecefc1ba6a5
md"""
## Q11

Suppose:

``A_{ij} \in \mathbb{Z}`` and

``\det(A) = \pm 1``

``A^{-1}_{ij} = \frac{C_{ji}}{\det(A)} \rightarrow |A^{-1}_{ij}| = C_{ji}``

Since:

``C_{ij} = (-1)^{i+j}\det(M_{ij})`` and we know that ``M_ij \in \mathbb{Z}``, ``A^{-1}_{ij} \in \mathbb{Z}``
"""

# ╔═╡ 8c0f2868-a834-42ba-a525-f18e0778e05e
q11ex = [
	1 2;
	3 5
]

# ╔═╡ 1fc159e6-4f03-4d2c-80b7-767452e50343
cofactor_determinant(q11ex)

# ╔═╡ 2815c57c-8bc0-46e9-90ca-ba851ea6ab06
q11ex\I

# ╔═╡ feacc1e3-d70d-4b64-bf2f-c388a8f0457a
md"""
## Q12

``\det(A)\det(A^{-1}) = \det(I) = \pm 1``
"""

# ╔═╡ 13010c70-a7ea-4156-98f3-5a48b1165b02
md"""
## Q13
"""

# ╔═╡ 6ccadc46-2798-443c-aa4d-57ec59f5d91d
q13A = [
	r 0 0;
	0 s 0;
	0 0 t
]

# ╔═╡ 10d4682e-2910-4b64-bd29-d19fa9fffca6
cofactor_determinant(q13A)

# ╔═╡ 9b5c5196-16cb-4dff-82d1-37e777855131
cramer(q13A)

# ╔═╡ bad2c32d-886e-4529-93d5-dfff34d223a8
md"""
## Q14
"""

# ╔═╡ 1ed39864-9adc-4f69-a7e5-bdd0aed46cee
q14L = [
	a 0 0;
	b c 0;
	d e f
]

# ╔═╡ 0b63e84a-d2d2-4aa4-97b6-331b18f8b3ab
all_cofactors(q14L)

# ╔═╡ 26ee7ecd-142c-4e06-8510-b28e5b14332a
cramer(q14L)

# ╔═╡ 5f621447-ec24-4b5c-8cc1-cfa942de300e
q14S = [
	a b d;
	b c e;
	d e f
]

# ╔═╡ 52f3244b-dedc-4b39-a2de-f9b7e0c95d1d
all_cofactors(q14S)

# ╔═╡ a2266a11-4c21-4af0-a6ab-731cabf0cfc9
cramer(q14S)

# ╔═╡ d8827bb9-e797-4e26-9141-163b6d355c5e
q14Q = [
	0 0 1;
	0 1 0;
	1 0 0
]

# ╔═╡ e50ed1f3-e16d-4545-ace1-aa5d766e38fc
all_cofactors(q14Q)

# ╔═╡ 01a3e0a6-6785-4a6b-8674-63c6b27a9943
md"""
``Q^{-1} = \frac{C^T}{\det(Q)}``

``Q^{-1} = \pm C^T``

``Q = \pm C``
"""

# ╔═╡ 3eb396ef-67ac-4abd-8cf7-27bf78d8d058
md"""
## Q15

Start with n = 2
4 cofactors containing 1 term each needing 2 multiplications

n = 5

25 cofactors containing 16 terms each needing 3! multiplications
"""

# ╔═╡ 909219e8-1fbe-490c-b0c0-f20a21221993
cofactor_multiplications(n) = n^2 * factorial(n-1) * 3

# ╔═╡ 2cc97818-d6a0-49c7-b01a-96ccaaca06ba
cofactor_multiplications(2)

# ╔═╡ 6b585b0e-ecf6-4420-976d-9d68cf994725
cofactor_multiplications(5)

# ╔═╡ ef7f21bc-3cd8-4ae5-8072-8e3e30d602b3
md"""
## Q16
"""

# ╔═╡ aae49b2e-17a8-40a7-9df5-053715836017
q16v = [3 2]

# ╔═╡ 34832b8e-224d-4cb0-8032-40aec5303088
q16w = [1 4]

# ╔═╡ 18bbe6e7-069a-483c-8859-1f17763811b4
q16A = vcat(q16v, q16w)

# ╔═╡ b905b0f3-0b76-48f6-bde7-e481a599e715
det(q16A)

# ╔═╡ 05518c71-1956-4e49-9159-b083bcb26234
# b
q16vw = q16w - q16v

# ╔═╡ a8f70b50-5180-4f6e-b6e2-535a2d1a2ac3
q16x = q16w + q16v

# ╔═╡ d03d4ff4-d589-4065-b4e4-4bcd67e90970
q16xw = q16x - q16v

# ╔═╡ 313a2270-9830-4d0b-aa80-53ca48095219
abs(0.5* det(vcat(q16vw, q16xw)))

# ╔═╡ 910d7a9c-4073-4f9e-b02f-fcc0d5dcfbc2
q16tri1 = Shape([0, q16v[1], q16w[1]], [0, q16v[2], q16w[2]])

# ╔═╡ 412573f6-4ae4-4fa9-abec-c85741949b18
plot(q16tri1, size = (1600, 1600))

# ╔═╡ 6e131d93-a4bc-4d76-8a69-82e077749c0e
q16y = q16vw - q16w

# ╔═╡ 842c13e4-6588-499e-83ca-8c0d32f4dfbc
abs(0.5* det(vcat(q16vw, q16y)))

# ╔═╡ 453f8ef9-82ca-4f41-960a-56d75fb64a97
q16tri2 = Shape([0, q16w[1], q16v[1]], [0, q16w[2], q16v[2]])

# ╔═╡ 476e32a2-ad79-4cb0-b337-53ee1fb7afa6
plot(q16tri1, size = (1600, 1600))

# ╔═╡ d061fb9f-7ac0-4791-bf49-792c40f64a33
md"""
## Q17
"""

# ╔═╡ 3072eccc-09d8-4c92-ae21-2efa6a99c445
q17u = [0 0 0] - [3 1 1]

# ╔═╡ 2a9dfc30-0034-499e-a287-2ece5d9ac9cd
q17v = [1 3 1]

# ╔═╡ f133fa5b-2678-475c-98c2-7c9af0fa114b
cross_product_length(q17u, q17v)

# ╔═╡ ed295efc-062e-4269-bdfd-645ea01e5dfc
cofactor_determinant([3 1 1; 1 3 1; 1 1 3])

# ╔═╡ 0874de30-df4b-4503-9994-6b5fe4809e34
md"""
## Q18
"""

# ╔═╡ c57972f9-320d-4d79-a463-bcf689d34044
q18u = [2 1]

# ╔═╡ 82cc004e-191b-4f35-80eb-e29ef79f5ab7
q18v = [3 4]

# ╔═╡ 600be562-68a2-4b01-9829-76cbad86663b
q18w = [0 5]

# ╔═╡ 1369b6be-925c-45ba-9da1-cbe4340ea61c
rat.(triangular_area(q18u, q18v, q18w))

# ╔═╡ 672da8be-cf1d-4d59-94f1-84a7d5102726
q18x = [-1 0]

# ╔═╡ 6588129e-802d-44ca-a2f6-2ccad132bba6
q18quad = Shape([x[1] for x in [q18u, q18v, q18w, q18x]], [x[2] for x in [q18u, q18v, q18w, q18x]])

# ╔═╡ b878c0c7-42ea-4198-a0d6-0aaf3f78afa8
plot(q18quad, size = (1600, 1600))

# ╔═╡ a1e9feec-fbf6-4b92-bd19-855a4bafb3ac
rat.(triangular_area(q18u, q18v, q18w)) + triangular_area(q18u, q18w, q18x)

# ╔═╡ 4aab9246-c59b-463f-bcc6-24a9a36e9d30
md"""
## Q19
"""

# ╔═╡ 23912767-7f91-46a3-922e-d818fd89ceef
q19u1 = [2 1]

# ╔═╡ e04efa2f-f2c7-4826-a2cf-d5cc97d7fd5b
q19v1 = [2 3]

# ╔═╡ 14d5c70b-4218-445f-afa2-fe95dddb76ef
triangular_area(q19u1, q19v1, q19u1 + q19v1) * 2

# ╔═╡ 65e46d42-842c-45ae-817e-5b9489a1a341
q19u2  = [2 2]

# ╔═╡ 652ebc69-fe28-4710-97da-02b23d4d588b
q19v2 = [1 3]

# ╔═╡ 623b53e0-dafd-4eeb-a0f0-9753cc9f2dda
triangular_area(q19u2, q19v2, q19u2 + q19v2) * 2

# ╔═╡ a8596118-a3d6-4b2d-a5e3-2f83d338627f
md"""
## Q20
"""

# ╔═╡ a57ffe79-4b4b-400a-9d70-1fa8e644621d
abs(cofactor_determinant(hadamard(4)))

# ╔═╡ 229f8c8a-752b-4180-a913-c8f4f7d50d31
md"""
## Q21

a) Largest possible value (because they are orthogonal to one another): ``V = \sqrt{L_1^2+ L_2^2 + L_3^2 + L_4^2}``
"""

# ╔═╡ add51acc-c66e-486f-92de-5d7ec67c4428
q21H = [
	1 1 1 1;
	1 1 1 -1;
	1 1 -1 -1;
	1 -1 -1 -1
]

# ╔═╡ 60cc5463-4b45-44b5-8928-7f52e30004ad
cofactor_determinant(q21H)

# ╔═╡ 6064ed69-edef-4c2d-89e2-ca42f0a5ba1c
cofactor_determinant(hadamard(4))

# ╔═╡ f06a8a06-dfb8-42b7-ad8d-057e88d846c4
md"""
## Q22

``x_1y_2 - x_2y_1`` is the same as the 2x2 determinant formula ``ad - bc``. Solving this geometrically would be difficult on the computer.
"""

# ╔═╡ afbadbdb-2134-4a48-b6a2-6695432f3ae0
md"""
## Q23

``\det(A) = ||a||||b||||c||`` since ``V = \det(A)``

``\det(A^TA) = ||a||^2||b||^2||c||^2`` Since the ``\det(A) = \det(A^T)`` and ``\det(A^TA) = \det(A^T)\det(A)``
"""

# ╔═╡ 60dcd6d7-6813-4f21-b1e9-211150ccd645
md"""
## Q24

``i`` is the height
"""

# ╔═╡ 0c17841f-6306-43cd-8b07-5f1cc9b7b262
q24u = [1 0 0]

# ╔═╡ d9659ea2-f314-4671-b934-6447e67746eb
q24v = [0 1 0]

# ╔═╡ b3788545-d090-4df9-9b16-75868f43e6c2
q24w = [2 3 4]

# ╔═╡ 52c6ec70-d84f-4c2c-aaa4-d5efb21ed1c4
q24A = vcat(q24u, q24v, q24w)

# ╔═╡ c850184f-beaa-4c9d-81fc-ee5360168433
cofactor_determinant(q24A) # This also the volume

# ╔═╡ 17c048a1-3f51-4bd2-bdaa-d57a282ac2fe
vector_length([2*i 3*j 4*k])

# ╔═╡ f6394348-9504-4186-bfb7-5873f37c6c10
cross_product_non_symbolic(q24u, q24v)

# ╔═╡ 383ebbc2-f3e6-4371-b4c9-743b17d5581e
triple_product(q24u, q24v, q24w)

# ╔═╡ 23aa5e54-2f80-442e-a6a8-4c2609266d8b
md"""
## Q25
"""

# ╔═╡ 6a1659b6-2270-40af-aa37-681bd7244acb
q25corners = [4, 8, 16]

# ╔═╡ 014c37c3-3754-439e-9aab-872349b146f5
md"""
The number of edges of a hypercube per corner is n, but need to divide by 2 because we  count twice.
"""

# ╔═╡ ceddaa60-811b-4d53-a5cd-9047b463be87
hypercube_edges(n) = n * 2 ^ (n - 1)

# ╔═╡ 1624de4f-cbc9-4001-bcff-2cd876a1b6dd
[hypercube_edges(i) for i=1:5]

# ╔═╡ aa97b07a-9c95-488d-b875-44f34145a296
md"""
We know 2-hypercube has 4 (n-1)-faces and, 3 has 6 (n - 1 faces). I suppose the number of cubes in a 4D-hypercube is 8. Coxeter shows the relation is

``E_{m,n} = 2^{n-m}\binom{n}{m}``
"""

# ╔═╡ 85e942d3-cfb5-46b1-b7ed-6eac55e6aa67
embedded_hypercubes(m,n) = Int(2^(n-m) * factorial(n) / (factorial(m) * factorial(n - m)))

# ╔═╡ a7467fdc-c534-472d-90b2-1f8a266b2a00
embedded_hypercubes(3, 4)

# ╔═╡ cd923378-8616-4238-820d-84adeefb26c9
md"""
The volume of a hypercube of ``2I`` is ``2^n``
"""

# ╔═╡ 0aec971d-04c6-4dee-bfbd-7109ad44fb4a
md"""
## Q26

``V_n = \int_0^1 V_{n-1}(x)dx``

``V_2 = \frac{x^2}{2}``

``V_3 = \int_0^1 V_2(x) = \frac{1}{2} \int_0^1 x^2 = \frac{1}{2} \left(\frac{1}{3} - 0 \right) = \frac{1}{6}``

``V_4 = \int_0^1 V_3(x) = \frac{1}{6} \int_0^1 x^3 = \frac{1}{6} \left(\frac{1}{4} - 0 \right) = \frac{1}{24}``

A 3-simplex has volume ``\frac{1}{6}``. A 4-simplex is volume ``\frac{1}{24}``
"""

# ╔═╡ 57fca92e-dea9-49fb-92de-f2f06897bb35
md"""
## Q27
"""

# ╔═╡ c4a3ea87-00d1-46fc-a086-d26101d60d12
q27J = [
	cos(θ) (-r * sin(θ));
	sin(θ) (r * cos(θ))
]

# ╔═╡ 208fc203-764a-4ae8-ba86-3706f0d3194a


# ╔═╡ 7738c8b7-965b-488c-a481-e16c5892e798
simplify(sqrt(sum(q27J[:,1] .^ 2))) # Same as 1

# ╔═╡ ad63cc9b-6157-412a-9907-0a219272b1f9
simplify(sqrt(sum(q27J[:,2] .^ 2))) # Same as r

# ╔═╡ b5efde02-3409-405f-bdf0-7fb701943a60
simplify(cofactor_determinant(q27J))

# ╔═╡ fd36b30f-0075-4260-8227-4f7ffc4ac635
md"""
## Q28
"""

# ╔═╡ eb6accdb-1bbf-4700-bd24-ff35a778be03
begin
	q28x = ρ * sin(ϕ) * cos(θ)
	q28y = ρ * sin(ϕ) * sin(θ)
	q28z = ρ * cos(ϕ)
end;

# ╔═╡ d728aa68-60f5-4577-af85-14cae44b0dcb
q28J = Symbolics.jacobian([q28x, q28y, q28z], [ρ, ϕ, θ], simplify=true)

# ╔═╡ b50f476d-c3db-4b2f-866e-a6fe3b7b01d6
simplify(cofactor_determinant(q28J), expand=true, thread_subtree_cutoff=1000)

# ╔═╡ 9c2ac996-bd54-4ecc-9ab1-c26c5b199143
md"""

``\det(J) = \sin^{3}\left( \phi \right) \cos^{2}\left( \theta \right) \rho^{2} + \sin^{3}\left( \phi \right) \sin^{2}\left( \theta \right) \rho^{2} + \cos^{2}\left( \phi \right) \cos^{2}\left( \theta \right) \rho^{2} \sin\left( \phi \right) + \cos^{2}\left( \phi \right) \sin^{2}\left( \theta \right) \rho^{2} \sin\left( \phi \right)``

``\det(J) = \rho^2\sin(\phi) (\sin^2(\phi)\cos^2(\theta) + \sin^2(\phi)\sin^2(\theta) + \cos^2(\phi)\cos^2(\theta) + \cos^2(\phi)\cos^2(\theta) )``

``\det(J) = \rho^2\sin(\phi) (\sin^2(\phi) + \cos^2(\phi))(\sin^2(\theta) + \cos^2(\theta)) = \rho^2\sin(\phi)``

"""

# ╔═╡ e66c2888-0dfe-4700-8b73-d601e5e1f319
Differential(q28x*q28y*q28z)(q28x)

# ╔═╡ 4837af43-641b-4b29-90cf-1ca2e320b90b
simplify(Differential(q28x*q28y*q28z)(ϕ))

# ╔═╡ b4e38353-7df1-4b27-877c-8f22942cbd42
md"""
## Q29
"""

# ╔═╡ 6a86d3ce-6529-4a81-98c2-e6c85bf96034
begin
	q29x = r * cos(θ)
	q29y = r * sin(θ)
end;

# ╔═╡ 34f929c0-a6c6-4594-9ef1-e427004cdf93
d29x = Differential(q29x)

# ╔═╡ 210b6890-f4cf-4688-bcc0-430d71b2d007
Symbolics.jacobian([q29x, q29y], [r, θ])

# ╔═╡ aac2023f-43ab-430f-9358-6602b8ee5563
begin
	q29r = x / cos(θ)
	q29θ = cos(x / r) ^ -1
end;

# ╔═╡ 28b12566-0f5c-4ed8-8193-a19ceac03bbe
q29Jinv = simplify.(two_by_two_inv(Symbolics.jacobian([q29x, q29y], [r, θ])), expand=true)

# ╔═╡ c6d689e9-b908-4d66-ad9c-33ef00f27e4d
simplify(det(q29Jinv)) # 1/r

# ╔═╡ e5b541fc-713c-4c81-9702-877dcfbc9a32
md"""
## Q30
"""

# ╔═╡ 732994bf-8e71-4041-aad3-0d791cd78498
begin 
	q30u = [0 0]
	q30v = [6 0]
	q30w = [1 4]
end;

# ╔═╡ 4765b8fd-a99e-4730-9425-b3875778c926
triangular_area(q30u, q30v, q30w) # Rotation shouldn't change area

# ╔═╡ 410fa9b7-aa75-4dbf-8f0e-447032b59dee
q30J = [
	cos(θ) -sin(θ);
	sin(θ) cos(θ)
]

# ╔═╡ a49bc510-2b84-441f-a060-007bad4c41aa
substitute.(q30J, (Dict([θ => 60 * π / 180]),))

# ╔═╡ 5bdbbaff-f9c4-445f-aa98-10684c48dd40
q30Jinv = simplify.(two_by_two_inv(q30J), expand=true)

# ╔═╡ 74e84aa9-f377-4eb3-8aca-b738514aa415
substitute.(q30Jinv, (Dict([θ => 60 * π / 180]),))

# ╔═╡ b7511dd3-cfc5-487c-a6b3-e56863682fef
det(substitute.(q30Jinv, (Dict([θ => 60 * π / 180]),)))

# ╔═╡ c2af082c-2f21-4b89-a568-c7e0880bcd94
md"""
## Q31
"""

# ╔═╡ d88fc0b6-f8cc-40b2-a53b-85e8c4870ae6
begin
	q31u = [2 4 0]
	q31v = [-1 3 0]
	q31w = [1 2 2]
end;

# ╔═╡ 9ba0b0cb-a133-4fc8-9e56-d8c56d552f2c
q31A = vector_length(cross_product_non_symbolic(q31u, q31v))

# ╔═╡ f8ef3e60-e48b-44b4-b039-0e3491b8ceb1
q31V = triple_product(q31u, q31v, q31w)

# ╔═╡ c058298c-2fad-42f1-828a-b17097388dd9
vector_length(q31w)

# ╔═╡ 2ab0da01-b912-4cc4-ae97-d2037efe3f76
acos(q31V / (q31A * vector_length(q31w)) * (π /180))

# ╔═╡ 530b13e1-62b1-48ff-a627-abb602f78567
q31V / q31A

# ╔═╡ ad50cf45-a168-49c8-aa03-d3a8a2995e67
md"""
## Q32
"""

# ╔═╡ 6c52a3de-0eaf-46f2-99d6-7afa653751c1
cofactor_determinant(vcat(q31u, q31v, q31w))

# ╔═╡ 3f6d8a71-6ac9-4748-ae3b-f550a54b35b8
md"""
## Q33
"""

# ╔═╡ e0f3eaa8-8ec8-4627-bb32-b0ec9f96992b
begin
	q33u = [u1 u2 u3]
	q33v = [v1 v2 v3]
	q33w = [w1 w2 w3]
end;

# ╔═╡ 4010aded-b562-4407-bd2f-aff0722fcac2
q33A = vcat(q33u, q33v, q33w)

# ╔═╡ 6778a88f-9455-4082-8122-8aa3ab39f102
cofactor_determinant(q33A)

# ╔═╡ cddbe5a2-d31a-47db-9ce2-7baa1b9b253f
md"""
## Q34

``(u \times v) \cdot w = \v{u_1 & u_2 & u_3 \\ v_1 & v_2 & v_3 \\ w_1 & w_2 & w_3} = \v{w_1 & w_2 & w_3 \\ u_1 & u_2 & u_3 \\ v_1 & v_2 & v_3}``


We need a triple product with an even number of permutations. We can start by looking for the triple product that places ``w`` on the top row:


``(w \times u) \cdot v``
"""

# ╔═╡ 3c1f401d-647d-40d1-911d-72b96f3b40ba
md"""
## Q35
"""

# ╔═╡ e78978da-e232-4739-860e-6e55b8d9a0ef
begin
	q35P = [1 0 -1]
	q35Q = [1 1 1]
	q35R = [2 2 1]
end;

# ╔═╡ bc868670-134c-4dd4-97a7-c5cee7e1994f
q35u = q35Q - q35P

# ╔═╡ be5e2e7c-cbef-4830-994f-556aaa43ad36
q35v = q35R - q35P

# ╔═╡ 36215464-d727-4d5e-b86b-c83c66d8b7b6
q35S = q35P + q35u + q35v

# ╔═╡ 17c43d95-b1fc-4578-886c-6b2848cca55b
vector_length(cross_product_non_symbolic(q35u, q35v))

# ╔═╡ 0c088b8f-24a0-4593-b3f7-29b8fb734455
q35O = [0 0 0]

# ╔═╡ f2648850-6266-4f1a-aa07-6687f031b972
q35w = q35O - q35P

# ╔═╡ 72e5c9d6-2a4b-47c5-be2a-f3a9680558b9
cross_product_non_symbolic(q35u, q35v)' * q35w'

# ╔═╡ 454a3771-d211-4e3f-950e-e600779fb9f0
q35T = q35P + q35w + q35u

# ╔═╡ db6ddeab-c433-4c87-a6f2-7903cb0d7eac
q35U = q35P + q35w + q35v

# ╔═╡ aeff9cba-b992-4b5b-857b-85f3c35f6cb5
q35V = q35P + q35u + q35v + q35w

# ╔═╡ 5e981696-4955-4951-8215-4f7eea3eb2c6
md"""
## Q36
"""

# ╔═╡ e8cba01b-a0b5-42ce-a661-2ca3dcef7821
begin
	q38u = [1 1 0]
	q38v = [1 2 1]
end;

# ╔═╡ 0a193abf-4051-474f-b9b1-94018746f78d
md"""
The deteriminant will be zero when ``\m{x & y & z} = au + bv``. 
"""

# ╔═╡ d901cace-3f8a-40f5-8deb-e4e0144257ab
substitute(cross_product(q38u, q38v), Dict(i => x, j => y, k => z))

# ╔═╡ b5831de7-fe61-4a10-bfbf-ca74676176c3
md"""
## Q37
"""

# ╔═╡ 9ac77af2-00cc-4848-8873-bc2df5db6dab
begin
	q39u = [2 3 1]
	q39v = [1 2 3]
end;

# ╔═╡ 5871e101-da65-4ee4-ba80-fb0bd4988e5d
substitute(cross_product(q39u, q39v), Dict(i => x, j => y, k => z))

# ╔═╡ e56e4e7f-15bd-45ea-b44a-51ec42ccbb84
md"""
## Q38

a) ``\det(2A) = 2^n\det(A)`` makes sense from a volumetric point of view if you consider that each row is an edge. Multiplying an edge by 2 doubles it's volumnes. Since you are using ``n`` edges, the doubling occurs ``n`` times.

b) ``\det(A) + \det(A) = \det(A + A)`` for a size 1 matrix.
"""

# ╔═╡ 7ea031f8-92c7-4683-88be-0bc91929cc79
det([3 1; 5 2] + [9 1; 4 2]) # counter-example for n=2

# ╔═╡ 8d1405e1-a945-401e-bccd-a98a87a2f344
det([3 1; 5 2])

# ╔═╡ a202bdb6-68d1-486b-8650-0fe54641cab7
det([9 1; 4 2])

# ╔═╡ 5b092c44-bfd0-4e4c-878a-e6663056e43e
md"""
## Q39

``A^{-1} = \frac{C^T}{\det(A)}``

``\det(A)A^{-1} = C^T``

``\det(\det(A)A^{-1}) = \det(C^T)``

``\det(A)^3  = \det(C^T)``

``\det(A) = \det(C^T)^{\frac{1}{3}}``

``A^{-1} = \frac{C^T}{\det(C^T)^{1/3}}``

``A = \frac{C^{-T}}{\det(C^T)^{1/3}}``
"""

# ╔═╡ 585734d1-a716-4531-ae9e-8a42a2199415
md"""
## Q40
"""

# ╔═╡ 19ec7991-308b-4268-bca8-4761453725e3
q40T = Tridiagonal(ones(4)*-1, ones(5) * 2, ones(4)*-1)

# ╔═╡ e564bf60-b461-4c92-b376-b50060199522
det(q40T)

# ╔═╡ 05aaef9b-6547-4ac4-9bd5-e68c46befbc0
det(q40T[1:2,1:2])

# ╔═╡ f0a599e7-c244-474f-916e-f273ea99b0fd
det(q40T[3:5,3:5])

# ╔═╡ 45842349-8c37-4769-990a-949f15962398
det(q40T[1:2,[1,3]])

# ╔═╡ 75101f0f-d568-4793-8b2b-f4a2e19d4484
det(q40T[3:5,[2,4,5]])

# ╔═╡ e9666346-6fa6-444d-840e-6351ca888da0
det(q40T[1:2,1:2])*det(q40T[3:5,3:5])-det(q40T[1:2,[1,3]])*det(q40T[3:5,[2,4,5]])

# ╔═╡ 67e7c046-9136-41c1-b558-3f3cc07453be
md"""
## Q41
"""

# ╔═╡ 7c169ecc-ca8b-468e-bf00-41b53f495513
q41A = [
	1 2 3;
	1 4 7
]

# ╔═╡ ffae66b7-126b-46c0-a48c-bd9374ee419c
q41B = q41A'

# ╔═╡ 6892be95-f9ed-4b5e-b63c-69f9e16c08ff
function cauchy_binet(A, B)
	(det(A[:,[1, 2]])*det(B[[1,2],:]) + det(A[:,[2,3]])*det(B[[2,3],:])) + det(A[:,[1,3]])*det(B[[1,3],:])
end

# ╔═╡ 29fc1e15-d4ad-43a6-ac0e-60c8552ecddf
cauchy_binet(q41A, q41B)

# ╔═╡ fa599fd3-31fa-47f2-b1c9-65da9d6d71a6
cofactor_determinant(q41A * q41B)

# ╔═╡ 7d939bbc-59f2-49ae-b17a-e9c6ec4cb46c
det(q41A[:,[1, 2]])

# ╔═╡ 1925090f-5c19-4f46-b4fd-6088a73ca85a
det(q41A[:,[1, 3]])

# ╔═╡ 715b5071-36d1-4c7d-9a52-ed2287457afa
det(q41A[:,[2, 3]])

# ╔═╡ a33e6f1c-6219-491a-84cf-7fa01ec08296
md"""
## Q42
"""

# ╔═╡ 06cb4f2c-37e9-4ee2-9118-0559343d5039
md"""
Consider a 5 by 5 tridiagonal matrix with cofactors `a_{11}` and `a_{12}`. `C_{11}` has the same number of values as the previous 4 by 4 matrix, which we know has 5 terms. `C_{12}` has two more zero values than ``C_{11}``. This causes the 4 x 4 matrix to degenerate into a 3 x 3 matrix, which has 3 values. So ``Terms(Tri(5)) = Terms(Tri(4)) + Terms(Tri(3))``


``C_{11} = Tri(4)``

``C_{12} = \m{x & 0 & 0 & 0 \\ 0 & x & y & 0 \\ 0 & z & x & y \\ 0 & 0 & z & x}``
"""

# ╔═╡ Cell order:
# ╠═26dfe76b-fbe3-4116-8129-a5538615df66
# ╠═20eb19f5-d2b1-425b-bd2e-4d142d28cd60
# ╠═923a3522-b6f8-434d-8bb9-85d8357df51d
# ╠═1644d02b-a035-4cc3-870e-5913b6bb3489
# ╠═763574f2-9d7c-441d-9fae-3d1253e39db7
# ╠═15495da5-b616-4540-8d60-6251ebcd8575
# ╠═39d59b60-bcc9-4c30-a08d-753913e68d91
# ╠═82ec58dd-0373-4438-86e4-4144bb0b0b34
# ╠═dc93074a-ebfe-4f4e-88f4-485f5c173d43
# ╠═43cc9465-3f44-4aad-a45d-5f1da3e67813
# ╠═58195e48-1e8a-4d16-bdc3-736ce4920a66
# ╠═fc2d2f95-edc8-41e6-8a10-b1173ddb51b0
# ╠═b207a078-2bfc-416b-afa2-6ae646647beb
# ╠═94599138-d357-43f1-876d-1073c0a49714
# ╠═5742bab2-f9fb-415b-a493-6e28b46d257b
# ╠═6a1227cb-3248-44ef-8ae1-d1842fb3d22f
# ╠═ab77210b-e692-48e6-a194-7b5562fa2e13
# ╠═4f1d0cb1-57a0-4dff-996c-372558de9aef
# ╠═fc5aef5a-1e0b-4316-86db-af2ab78110dc
# ╠═98422293-89b8-4f80-830b-9c360f52c42d
# ╠═c2b022d8-2ecf-43cf-844d-4b9e4bc198c4
# ╠═c0b1135e-90b3-4308-b668-4ddd3c5b97fe
# ╠═c6ff16a6-96b8-4f45-af2d-05aca7d97d84
# ╠═3272d2ca-50c0-4c2f-a363-935af2bf00b3
# ╠═52a08df7-5625-4833-9ab7-793534697114
# ╠═70bd870f-cc32-4840-8ea1-f54bb02855e4
# ╠═56fc6a16-5502-4116-a507-1559b6672786
# ╠═1b3b179c-2aa2-4fa4-9eaa-0ccfcd64ed4b
# ╠═e2993de4-bd70-11eb-0ad0-75d1f93d8a80
# ╠═0f017eb0-a10e-420a-9d68-4a77348b1c26
# ╠═802f6c24-5657-4e27-bfab-d340ee6111ff
# ╠═c41c25c3-3b62-43e1-aada-b90c3de9efed
# ╠═9e6f0924-ce5e-4827-a4f9-d8f47d2f2fea
# ╠═c064ad99-9ce2-4981-8607-04c1c7a5d85a
# ╠═ca6821e1-bff3-4aa6-b9c7-922f258ae55c
# ╠═9660772f-bbe6-4e63-90d5-212862840c2b
# ╠═f3bc7f53-00bb-4ede-b353-189fd41c60fd
# ╠═b54a70e9-a45f-4156-b51b-b8fb4e030c32
# ╠═d230dc4e-6ef2-4919-95ed-9f2603bcfd9f
# ╠═2e725c38-318c-479b-901e-98733856df60
# ╠═ecdcac60-2410-4038-af98-a51d67f3e266
# ╠═1cf8bc16-a088-4f82-b414-5c381376e744
# ╠═56b0a487-a01a-42dd-9382-db38c5d4e001
# ╠═e2344c48-5145-45d3-a26b-7128e20f9316
# ╠═4ff64552-e76d-4d8f-89d1-d8b3d0f98aae
# ╠═da6ed2e3-0a1a-4126-81a8-2f077823ce87
# ╠═1b2b7814-6550-454c-b323-8b776991940a
# ╠═3e2c4dc5-87b4-479d-9a79-f7e986cdf412
# ╠═ffda5ac4-0964-45f0-b0ed-45b99c5cc883
# ╠═1ead8e56-5791-4951-ad89-d26ee94e3bec
# ╠═c4f284af-c4ee-4e3b-84ac-e2d01d2a2f58
# ╠═7f5612c7-1e5a-4206-8f0b-927b78ba6625
# ╠═9f3a3a10-d447-49ae-9ed2-22a0d96235a5
# ╠═96387e53-dc09-403a-b900-7604ff6ecbe3
# ╠═18fa8409-ac1c-4d56-b6d9-26f0ad30214b
# ╠═f0044b20-a546-4c88-8f0b-6a16c98ac33a
# ╠═48e7e5d3-305b-4591-9cab-d8b70ee9b595
# ╠═483c1ff0-85f2-469e-ba15-80910ba80fd4
# ╠═8861908e-0442-426b-940f-3b4ba9996a18
# ╠═bcd49c46-388a-4d64-846e-1f81d38c2f6a
# ╠═6d661775-d16d-4319-bfc7-c905b5fcc242
# ╠═d72684f0-7dbf-4875-b974-f41b6f530b2f
# ╠═0e64596d-8408-4d36-b7c8-41bcca851bb9
# ╠═c6727867-89a4-4cb6-9858-eb4acea7b338
# ╠═ead0878d-fcea-4666-ba11-b665e7a49880
# ╠═29fc5286-003d-41f1-9b68-404ff0537c40
# ╠═fbddb9b0-ac15-4fee-bb18-c15d2951a083
# ╠═7bbf2f1d-c963-4f9d-ae88-ea0fa80f10e8
# ╠═dc319874-3405-4a50-a7b6-a93dd9b3ca0a
# ╠═962a0d8a-bff8-49ac-bf56-5aec93f154ff
# ╠═5835b954-24d0-4d2c-84af-e020796c5e5f
# ╠═52efbead-2202-4fd3-b3da-1fa9580fbe16
# ╠═66041ea1-6876-47f8-b3db-42a722b4531b
# ╠═30939bc8-6eca-48bd-8e11-bb68de9428fc
# ╠═a2a52c1e-ec88-488e-b1c4-185f15a0bea5
# ╠═19f93fb4-6b2a-4c19-b755-6ecefc1ba6a5
# ╠═8c0f2868-a834-42ba-a525-f18e0778e05e
# ╠═1fc159e6-4f03-4d2c-80b7-767452e50343
# ╠═2815c57c-8bc0-46e9-90ca-ba851ea6ab06
# ╠═feacc1e3-d70d-4b64-bf2f-c388a8f0457a
# ╠═13010c70-a7ea-4156-98f3-5a48b1165b02
# ╠═6ccadc46-2798-443c-aa4d-57ec59f5d91d
# ╠═10d4682e-2910-4b64-bd29-d19fa9fffca6
# ╠═9b5c5196-16cb-4dff-82d1-37e777855131
# ╠═bad2c32d-886e-4529-93d5-dfff34d223a8
# ╠═1ed39864-9adc-4f69-a7e5-bdd0aed46cee
# ╠═0b63e84a-d2d2-4aa4-97b6-331b18f8b3ab
# ╠═26ee7ecd-142c-4e06-8510-b28e5b14332a
# ╠═5f621447-ec24-4b5c-8cc1-cfa942de300e
# ╠═52f3244b-dedc-4b39-a2de-f9b7e0c95d1d
# ╠═a2266a11-4c21-4af0-a6ab-731cabf0cfc9
# ╠═d8827bb9-e797-4e26-9141-163b6d355c5e
# ╠═e50ed1f3-e16d-4545-ace1-aa5d766e38fc
# ╠═01a3e0a6-6785-4a6b-8674-63c6b27a9943
# ╠═3eb396ef-67ac-4abd-8cf7-27bf78d8d058
# ╠═909219e8-1fbe-490c-b0c0-f20a21221993
# ╠═2cc97818-d6a0-49c7-b01a-96ccaaca06ba
# ╠═6b585b0e-ecf6-4420-976d-9d68cf994725
# ╠═ef7f21bc-3cd8-4ae5-8072-8e3e30d602b3
# ╠═aae49b2e-17a8-40a7-9df5-053715836017
# ╠═34832b8e-224d-4cb0-8032-40aec5303088
# ╠═18bbe6e7-069a-483c-8859-1f17763811b4
# ╠═b905b0f3-0b76-48f6-bde7-e481a599e715
# ╠═05518c71-1956-4e49-9159-b083bcb26234
# ╠═a8f70b50-5180-4f6e-b6e2-535a2d1a2ac3
# ╠═d03d4ff4-d589-4065-b4e4-4bcd67e90970
# ╠═313a2270-9830-4d0b-aa80-53ca48095219
# ╠═910d7a9c-4073-4f9e-b02f-fcc0d5dcfbc2
# ╠═412573f6-4ae4-4fa9-abec-c85741949b18
# ╠═6e131d93-a4bc-4d76-8a69-82e077749c0e
# ╠═842c13e4-6588-499e-83ca-8c0d32f4dfbc
# ╠═453f8ef9-82ca-4f41-960a-56d75fb64a97
# ╠═476e32a2-ad79-4cb0-b337-53ee1fb7afa6
# ╠═d061fb9f-7ac0-4791-bf49-792c40f64a33
# ╠═3072eccc-09d8-4c92-ae21-2efa6a99c445
# ╠═2a9dfc30-0034-499e-a287-2ece5d9ac9cd
# ╠═f133fa5b-2678-475c-98c2-7c9af0fa114b
# ╠═ed295efc-062e-4269-bdfd-645ea01e5dfc
# ╠═0874de30-df4b-4503-9994-6b5fe4809e34
# ╠═c57972f9-320d-4d79-a463-bcf689d34044
# ╠═82cc004e-191b-4f35-80eb-e29ef79f5ab7
# ╠═600be562-68a2-4b01-9829-76cbad86663b
# ╠═1369b6be-925c-45ba-9da1-cbe4340ea61c
# ╠═672da8be-cf1d-4d59-94f1-84a7d5102726
# ╠═b878c0c7-42ea-4198-a0d6-0aaf3f78afa8
# ╠═6588129e-802d-44ca-a2f6-2ccad132bba6
# ╠═a1e9feec-fbf6-4b92-bd19-855a4bafb3ac
# ╠═4aab9246-c59b-463f-bcc6-24a9a36e9d30
# ╠═23912767-7f91-46a3-922e-d818fd89ceef
# ╠═e04efa2f-f2c7-4826-a2cf-d5cc97d7fd5b
# ╠═14d5c70b-4218-445f-afa2-fe95dddb76ef
# ╠═65e46d42-842c-45ae-817e-5b9489a1a341
# ╠═652ebc69-fe28-4710-97da-02b23d4d588b
# ╠═623b53e0-dafd-4eeb-a0f0-9753cc9f2dda
# ╠═a8596118-a3d6-4b2d-a5e3-2f83d338627f
# ╠═a57ffe79-4b4b-400a-9d70-1fa8e644621d
# ╠═229f8c8a-752b-4180-a913-c8f4f7d50d31
# ╠═add51acc-c66e-486f-92de-5d7ec67c4428
# ╠═60cc5463-4b45-44b5-8928-7f52e30004ad
# ╠═6064ed69-edef-4c2d-89e2-ca42f0a5ba1c
# ╠═f06a8a06-dfb8-42b7-ad8d-057e88d846c4
# ╠═afbadbdb-2134-4a48-b6a2-6695432f3ae0
# ╠═60dcd6d7-6813-4f21-b1e9-211150ccd645
# ╠═0c17841f-6306-43cd-8b07-5f1cc9b7b262
# ╠═d9659ea2-f314-4671-b934-6447e67746eb
# ╠═b3788545-d090-4df9-9b16-75868f43e6c2
# ╠═52c6ec70-d84f-4c2c-aaa4-d5efb21ed1c4
# ╠═c850184f-beaa-4c9d-81fc-ee5360168433
# ╠═17c048a1-3f51-4bd2-bdaa-d57a282ac2fe
# ╠═f6394348-9504-4186-bfb7-5873f37c6c10
# ╠═383ebbc2-f3e6-4371-b4c9-743b17d5581e
# ╠═23aa5e54-2f80-442e-a6a8-4c2609266d8b
# ╠═6a1659b6-2270-40af-aa37-681bd7244acb
# ╠═014c37c3-3754-439e-9aab-872349b146f5
# ╠═ceddaa60-811b-4d53-a5cd-9047b463be87
# ╠═1624de4f-cbc9-4001-bcff-2cd876a1b6dd
# ╠═aa97b07a-9c95-488d-b875-44f34145a296
# ╠═85e942d3-cfb5-46b1-b7ed-6eac55e6aa67
# ╠═a7467fdc-c534-472d-90b2-1f8a266b2a00
# ╠═cd923378-8616-4238-820d-84adeefb26c9
# ╠═0aec971d-04c6-4dee-bfbd-7109ad44fb4a
# ╠═57fca92e-dea9-49fb-92de-f2f06897bb35
# ╠═c4a3ea87-00d1-46fc-a086-d26101d60d12
# ╠═208fc203-764a-4ae8-ba86-3706f0d3194a
# ╠═7738c8b7-965b-488c-a481-e16c5892e798
# ╠═ad63cc9b-6157-412a-9907-0a219272b1f9
# ╠═b5efde02-3409-405f-bdf0-7fb701943a60
# ╠═fd36b30f-0075-4260-8227-4f7ffc4ac635
# ╠═eb6accdb-1bbf-4700-bd24-ff35a778be03
# ╠═d728aa68-60f5-4577-af85-14cae44b0dcb
# ╠═b50f476d-c3db-4b2f-866e-a6fe3b7b01d6
# ╠═9c2ac996-bd54-4ecc-9ab1-c26c5b199143
# ╠═e66c2888-0dfe-4700-8b73-d601e5e1f319
# ╠═4837af43-641b-4b29-90cf-1ca2e320b90b
# ╠═b4e38353-7df1-4b27-877c-8f22942cbd42
# ╠═6a86d3ce-6529-4a81-98c2-e6c85bf96034
# ╠═34f929c0-a6c6-4594-9ef1-e427004cdf93
# ╠═210b6890-f4cf-4688-bcc0-430d71b2d007
# ╠═aac2023f-43ab-430f-9358-6602b8ee5563
# ╠═28b12566-0f5c-4ed8-8193-a19ceac03bbe
# ╠═c6d689e9-b908-4d66-ad9c-33ef00f27e4d
# ╠═e5b541fc-713c-4c81-9702-877dcfbc9a32
# ╠═732994bf-8e71-4041-aad3-0d791cd78498
# ╠═4765b8fd-a99e-4730-9425-b3875778c926
# ╠═410fa9b7-aa75-4dbf-8f0e-447032b59dee
# ╠═a49bc510-2b84-441f-a060-007bad4c41aa
# ╠═5bdbbaff-f9c4-445f-aa98-10684c48dd40
# ╠═74e84aa9-f377-4eb3-8aca-b738514aa415
# ╠═b7511dd3-cfc5-487c-a6b3-e56863682fef
# ╠═c2af082c-2f21-4b89-a568-c7e0880bcd94
# ╠═d88fc0b6-f8cc-40b2-a53b-85e8c4870ae6
# ╠═9ba0b0cb-a133-4fc8-9e56-d8c56d552f2c
# ╠═f8ef3e60-e48b-44b4-b039-0e3491b8ceb1
# ╠═c058298c-2fad-42f1-828a-b17097388dd9
# ╠═2ab0da01-b912-4cc4-ae97-d2037efe3f76
# ╠═530b13e1-62b1-48ff-a627-abb602f78567
# ╠═ad50cf45-a168-49c8-aa03-d3a8a2995e67
# ╠═6c52a3de-0eaf-46f2-99d6-7afa653751c1
# ╠═3f6d8a71-6ac9-4748-ae3b-f550a54b35b8
# ╠═e0f3eaa8-8ec8-4627-bb32-b0ec9f96992b
# ╠═4010aded-b562-4407-bd2f-aff0722fcac2
# ╠═6778a88f-9455-4082-8122-8aa3ab39f102
# ╠═cddbe5a2-d31a-47db-9ce2-7baa1b9b253f
# ╠═3c1f401d-647d-40d1-911d-72b96f3b40ba
# ╠═e78978da-e232-4739-860e-6e55b8d9a0ef
# ╠═bc868670-134c-4dd4-97a7-c5cee7e1994f
# ╠═be5e2e7c-cbef-4830-994f-556aaa43ad36
# ╠═36215464-d727-4d5e-b86b-c83c66d8b7b6
# ╠═17c43d95-b1fc-4578-886c-6b2848cca55b
# ╠═0c088b8f-24a0-4593-b3f7-29b8fb734455
# ╠═f2648850-6266-4f1a-aa07-6687f031b972
# ╠═72e5c9d6-2a4b-47c5-be2a-f3a9680558b9
# ╠═454a3771-d211-4e3f-950e-e600779fb9f0
# ╠═db6ddeab-c433-4c87-a6f2-7903cb0d7eac
# ╠═aeff9cba-b992-4b5b-857b-85f3c35f6cb5
# ╠═5e981696-4955-4951-8215-4f7eea3eb2c6
# ╠═e8cba01b-a0b5-42ce-a661-2ca3dcef7821
# ╠═0a193abf-4051-474f-b9b1-94018746f78d
# ╠═d901cace-3f8a-40f5-8deb-e4e0144257ab
# ╠═b5831de7-fe61-4a10-bfbf-ca74676176c3
# ╠═9ac77af2-00cc-4848-8873-bc2df5db6dab
# ╠═5871e101-da65-4ee4-ba80-fb0bd4988e5d
# ╠═e56e4e7f-15bd-45ea-b44a-51ec42ccbb84
# ╠═7ea031f8-92c7-4683-88be-0bc91929cc79
# ╠═8d1405e1-a945-401e-bccd-a98a87a2f344
# ╠═a202bdb6-68d1-486b-8650-0fe54641cab7
# ╠═5b092c44-bfd0-4e4c-878a-e6663056e43e
# ╠═585734d1-a716-4531-ae9e-8a42a2199415
# ╠═19ec7991-308b-4268-bca8-4761453725e3
# ╠═e564bf60-b461-4c92-b376-b50060199522
# ╠═05aaef9b-6547-4ac4-9bd5-e68c46befbc0
# ╠═f0a599e7-c244-474f-916e-f273ea99b0fd
# ╠═45842349-8c37-4769-990a-949f15962398
# ╠═75101f0f-d568-4793-8b2b-f4a2e19d4484
# ╠═e9666346-6fa6-444d-840e-6351ca888da0
# ╠═67e7c046-9136-41c1-b558-3f3cc07453be
# ╠═7c169ecc-ca8b-468e-bf00-41b53f495513
# ╠═ffae66b7-126b-46c0-a48c-bd9374ee419c
# ╠═6892be95-f9ed-4b5e-b63c-69f9e16c08ff
# ╠═29fc1e15-d4ad-43a6-ac0e-60c8552ecddf
# ╠═fa599fd3-31fa-47f2-b1c9-65da9d6d71a6
# ╠═7d939bbc-59f2-49ae-b17a-e9c6ec4cb46c
# ╠═1925090f-5c19-4f46-b4fd-6088a73ca85a
# ╠═715b5071-36d1-4c7d-9a52-ed2287457afa
# ╠═a33e6f1c-6219-491a-84cf-7fa01ec08296
# ╠═06cb4f2c-37e9-4ee2-9118-0559343d5039
