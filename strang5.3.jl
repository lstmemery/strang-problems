### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 26dfe76b-fbe3-4116-8129-a5538615df66
using LinearAlgebra

# ╔═╡ 20eb19f5-d2b1-425b-bd2e-4d142d28cd60
using LaTeXStrings

# ╔═╡ 923a3522-b6f8-434d-8bb9-85d8357df51d
using Symbolics

# ╔═╡ 15495da5-b616-4540-8d60-6251ebcd8575
using Plots

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
	@variables x a b c d e f g h i j k l m n o p q r s t u v w y z D
end;

# ╔═╡ 98422293-89b8-4f80-830b-9c360f52c42d
function cross_product(u, v)
	A = vcat([i j k], u,v)
	cofactor_determinant(A)
end

# ╔═╡ 3272d2ca-50c0-4c2f-a363-935af2bf00b3
function cross_product_length(u,v)
	cross = cross_product(u,v)
	sqrt(sum(values(Symbolics.value(cross).dict) .^ 2))
end

# ╔═╡ 52a08df7-5625-4833-9ab7-793534697114
function triangular_area(u, v, w)
	tri_matrix = hvcat((3, 3, 3), u..., 1, v..., 1, w..., 1)
	cofactor_determinant(tri_matrix) / 2
end

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
q19u = [2 1]

# ╔═╡ e04efa2f-f2c7-4826-a2cf-d5cc97d7fd5b
q19v = [2 3]

# ╔═╡ a8596118-a3d6-4b2d-a5e3-2f83d338627f
md"""
## Q20
"""

# ╔═╡ 229f8c8a-752b-4180-a913-c8f4f7d50d31
md"""
## Q21
"""

# ╔═╡ f06a8a06-dfb8-42b7-ad8d-057e88d846c4
md"""
## Q22
"""

# ╔═╡ afbadbdb-2134-4a48-b6a2-6695432f3ae0
md"""
## Q23
"""

# ╔═╡ 60dcd6d7-6813-4f21-b1e9-211150ccd645
md"""
## Q24
"""

# ╔═╡ 23aa5e54-2f80-442e-a6a8-4c2609266d8b
md"""
## Q25
"""

# ╔═╡ 0aec971d-04c6-4dee-bfbd-7109ad44fb4a
md"""
## Q26
"""

# ╔═╡ 57fca92e-dea9-49fb-92de-f2f06897bb35
md"""
## Q27
"""

# ╔═╡ fd36b30f-0075-4260-8227-4f7ffc4ac635
md"""
## Q28
"""

# ╔═╡ b4e38353-7df1-4b27-877c-8f22942cbd42
md"""
## Q29
"""

# ╔═╡ e5b541fc-713c-4c81-9702-877dcfbc9a32
md"""
## Q30
"""

# ╔═╡ c2af082c-2f21-4b89-a568-c7e0880bcd94
md"""
## Q31
"""

# ╔═╡ ad50cf45-a168-49c8-aa03-d3a8a2995e67
md"""
## Q32
"""

# ╔═╡ 3f6d8a71-6ac9-4748-ae3b-f550a54b35b8
md"""
## Q33
"""

# ╔═╡ cddbe5a2-d31a-47db-9ce2-7baa1b9b253f
md"""
## Q34
"""

# ╔═╡ 3c1f401d-647d-40d1-911d-72b96f3b40ba
md"""
## Q35
"""

# ╔═╡ 5e981696-4955-4951-8215-4f7eea3eb2c6
md"""
## Q36
"""

# ╔═╡ b5831de7-fe61-4a10-bfbf-ca74676176c3
md"""
## Q37
"""

# ╔═╡ e56e4e7f-15bd-45ea-b44a-51ec42ccbb84
md"""
## Q38
"""

# ╔═╡ 5b092c44-bfd0-4e4c-878a-e6663056e43e
md"""
## Q39
"""

# ╔═╡ 585734d1-a716-4531-ae9e-8a42a2199415
md"""
## Q40
"""

# ╔═╡ 67e7c046-9136-41c1-b558-3f3cc07453be
md"""
## Q41
"""

# ╔═╡ a33e6f1c-6219-491a-84cf-7fa01ec08296
md"""
## Q42
"""

# ╔═╡ 06cb4f2c-37e9-4ee2-9118-0559343d5039


# ╔═╡ Cell order:
# ╠═26dfe76b-fbe3-4116-8129-a5538615df66
# ╠═20eb19f5-d2b1-425b-bd2e-4d142d28cd60
# ╠═923a3522-b6f8-434d-8bb9-85d8357df51d
# ╠═15495da5-b616-4540-8d60-6251ebcd8575
# ╠═82ec58dd-0373-4438-86e4-4144bb0b0b34
# ╠═dc93074a-ebfe-4f4e-88f4-485f5c173d43
# ╠═43cc9465-3f44-4aad-a45d-5f1da3e67813
# ╠═58195e48-1e8a-4d16-bdc3-736ce4920a66
# ╠═fc2d2f95-edc8-41e6-8a10-b1173ddb51b0
# ╠═b207a078-2bfc-416b-afa2-6ae646647beb
# ╠═5742bab2-f9fb-415b-a493-6e28b46d257b
# ╠═6a1227cb-3248-44ef-8ae1-d1842fb3d22f
# ╠═ab77210b-e692-48e6-a194-7b5562fa2e13
# ╠═4f1d0cb1-57a0-4dff-996c-372558de9aef
# ╠═fc5aef5a-1e0b-4316-86db-af2ab78110dc
# ╠═98422293-89b8-4f80-830b-9c360f52c42d
# ╠═3272d2ca-50c0-4c2f-a363-935af2bf00b3
# ╠═52a08df7-5625-4833-9ab7-793534697114
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
# ╠═a8596118-a3d6-4b2d-a5e3-2f83d338627f
# ╠═229f8c8a-752b-4180-a913-c8f4f7d50d31
# ╠═f06a8a06-dfb8-42b7-ad8d-057e88d846c4
# ╠═afbadbdb-2134-4a48-b6a2-6695432f3ae0
# ╠═60dcd6d7-6813-4f21-b1e9-211150ccd645
# ╠═23aa5e54-2f80-442e-a6a8-4c2609266d8b
# ╠═0aec971d-04c6-4dee-bfbd-7109ad44fb4a
# ╠═57fca92e-dea9-49fb-92de-f2f06897bb35
# ╠═fd36b30f-0075-4260-8227-4f7ffc4ac635
# ╠═b4e38353-7df1-4b27-877c-8f22942cbd42
# ╠═e5b541fc-713c-4c81-9702-877dcfbc9a32
# ╠═c2af082c-2f21-4b89-a568-c7e0880bcd94
# ╠═ad50cf45-a168-49c8-aa03-d3a8a2995e67
# ╠═3f6d8a71-6ac9-4748-ae3b-f550a54b35b8
# ╠═cddbe5a2-d31a-47db-9ce2-7baa1b9b253f
# ╠═3c1f401d-647d-40d1-911d-72b96f3b40ba
# ╠═5e981696-4955-4951-8215-4f7eea3eb2c6
# ╠═b5831de7-fe61-4a10-bfbf-ca74676176c3
# ╠═e56e4e7f-15bd-45ea-b44a-51ec42ccbb84
# ╠═5b092c44-bfd0-4e4c-878a-e6663056e43e
# ╠═585734d1-a716-4531-ae9e-8a42a2199415
# ╠═67e7c046-9136-41c1-b558-3f3cc07453be
# ╠═a33e6f1c-6219-491a-84cf-7fa01ec08296
# ╠═06cb4f2c-37e9-4ee2-9118-0559343d5039
