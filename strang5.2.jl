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

# ╔═╡ 111fa6a0-3d84-420e-b546-59075a03a5c2
@variables x a b c d

# ╔═╡ 5b4a23b6-17ee-4ce9-814f-cb8dd2489a58
rat(a) = rationalize(a, tol=1e-4)

# ╔═╡ 49f2380f-af3b-4b20-afde-1b793499cdec
md"""
``\newcommand\m[1]{\begin{bmatrix}#1\end{bmatrix}}``
"""

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

# ╔═╡ 33714c08-90c6-4044-ad80-c5649024f47f
q1B = [
	1 2 3;
	4 4 4;
	5 6 7
]

# ╔═╡ a09e4c71-a672-41f8-8989-14f583da0006
q1C = [
	1 1 1;
	1 1 0;
	1 0 0
]

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

# ╔═╡ 660eea34-7f99-4214-bccd-1ee0a299405a
q2B = [
	1 2 3;
	4 5 6;
	7 8 9
]

# ╔═╡ 61948f65-5a96-415d-88f9-d05169b7d775
q2C = vcat(hcat(q2A,zeros(3, 3)), hcat(zeros(3,3),q2A))

# ╔═╡ 478ce283-3b5c-4a19-8aac-a8e7678928b0
q2D = vcat(hcat(q2A,zeros(3, 3)), hcat(zeros(3,3),q2B))

# ╔═╡ c4b032d4-9ab4-4ebd-baed-8fe2ca0a9a64
md"""
## Q3
"""

# ╔═╡ a041250a-6a9a-47d8-9dbe-ee562722a53f
q3A = [
	x x x;
	0 0 x;
	0 0 x
]

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

# ╔═╡ 3d512797-1045-4c78-9a17-c6aa0d7a80a9
q4B = [
	1 0 0 2;
	0 3 4 5;
	5 4 0 3;
	2 0 0 1
]

# ╔═╡ 672ce5a4-0f6a-48ce-bb54-dad1371ea60c
md"""
## Q5
"""

# ╔═╡ 03b9475f-efdf-4897-8797-e8498f343062
md"""
## Q6
"""

# ╔═╡ d2f3b244-dc7f-4a51-99bb-cfeffb64d21b
md"""
## Q7
"""

# ╔═╡ 9361ae3c-ab6a-4c9b-8eee-bac81c78eb98
md"""
## Q8
"""

# ╔═╡ c3bfaea3-db75-456a-ae17-9bdfc5cacdef
md"""
## Q9
"""

# ╔═╡ 6a133cee-048e-4c26-9c11-be71a22edddc
md"""
## Q10
"""

# ╔═╡ fbe30313-f9e9-4151-b426-df61b56dc28f
md"""
## Q11
"""

# ╔═╡ 433b3e27-fd82-4b19-ad90-102618fcc381
q11A = [
	a b;
	c d
]

# ╔═╡ eddef1e3-ffff-400e-af18-88cfe67a81b2
q11B = [
	1 2 3;
	4 5 6;
	7 0 0
]

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

# ╔═╡ 58bd585b-ac46-40d2-b51f-10c0f56225d2
q13C3 = tridiagonal13(3)

# ╔═╡ e4570afd-c1b1-49e5-a4a8-29a092a45c61
q13C4 = tridiagonal13(4)

# ╔═╡ b4dde61d-44fc-42ed-9e33-4efcbfcc9bc3
md"""
## Q14
"""

# ╔═╡ 221f3f2d-442e-48bd-aa2c-5fba5e069844
md"""
## Q15
"""

# ╔═╡ ac3a60fc-4ac4-40a7-9c68-a0bf860dea01
q15E1 = [1]

# ╔═╡ 44fff8ac-9940-4cb4-9ab2-ee5ea9ced918
tridiagonal15(size) = Tridiagonal(ones(size-1), ones(size), ones(size-1))

# ╔═╡ 2fa48fd0-9047-4308-b31d-f0ea5dc4b91a
q15E2 = tridiagonal15(2)

# ╔═╡ 2e6cdf1e-295b-4554-8c13-9737eabbb59a
q15E3 = tridiagonal15(3)

# ╔═╡ f9710878-110e-4ef7-a3a1-ed47f22d7dcc
q15E4 = tridiagonal15(4)

# ╔═╡ c5a41df2-b049-4ec7-947d-22c97ef65b38
md"""
## Q16
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
"""

# ╔═╡ 292f3d65-7669-4326-9597-73313345c3e9
tridiagonal17(size) = Tridiagonal(-1*ones(size-1), vec(hcat(1, ones(size-1)*2...)), -1*ones(size-1))

# ╔═╡ 25e06f36-f583-41c1-8aca-a92bcc780dbf
q17B4 = tridiagonal17(4)

# ╔═╡ 64d53601-3b03-4155-b609-b3864fcc5b31
q17B3 = tridiagonal17(3)

# ╔═╡ 14ecaf00-2497-45b6-9c38-65644ff59745
q17B2 = tridiagonal17(2)

# ╔═╡ 07105032-916f-4caa-850f-fdacd8a1b971
md"""
## Q18
"""

# ╔═╡ c3428cf4-cf48-4b9d-ba3c-da8187018bdd
md"""
## Q19
"""

# ╔═╡ e2c7912a-d980-466a-8843-ae13893839b6
q19V4 = [
	1 a a^2 a^3;
	1 b b^2 b^3;
	1 c c^2 c^3;
	a x x^2 x^3
]

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

# ╔═╡ 84449c10-438a-44b8-ac1e-dd617296dc8a
md"""
## Q22
"""

# ╔═╡ e4eb0ee4-76eb-4b63-ae81-6efaf90c5338
constructq22(size) = Tridiagonal(ones(size-1), vec(hcat(2, ones(size-1)*3...)), ones(size-1))

# ╔═╡ 36dd5725-1589-4f29-ae6c-9f4ebe223e83
q22S3 = constructq22(3)

# ╔═╡ a886d29e-bce2-4316-860d-3063df707f36
md"""
## Q23
"""

# ╔═╡ 4447861d-56f6-44ae-91d6-87a393be0a43
md"""
## Q24
"""

# ╔═╡ eb0c5ef2-28d5-431f-8d22-7db0564c24a3
md"""
## Q25
"""

# ╔═╡ 6a6ec27c-4cd1-429b-b881-f92552550118
md"""
## Q26
"""

# ╔═╡ 63e2463d-fa27-44cf-a2cc-fc53caa019f6
md"""
## Q27
"""

# ╔═╡ 13b14998-2bf9-4e5d-b902-d480a5d9711c
md"""
## Q28
"""

# ╔═╡ 24ddccd5-8787-4c58-8c62-4dca6787762f
md"""
## Q29
"""

# ╔═╡ 03399151-21e9-4e5f-aa6e-a7813ca47a15
md"""
## Q30
"""

# ╔═╡ cb514a8a-123c-4839-a17d-8d1fb41430ab
constructq30(size) = Tridiagonal(ones(size-1)*-1, ones(size-1) * 2, ones(size-1) * -1)

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

# ╔═╡ 8baed9ad-8824-41f8-ab8a-c2acf60edb46
md"""
## Q32
"""

# ╔═╡ 0a8fad3e-759e-4b8b-a3a0-466574feb9c9
md"""
## Q33
"""

# ╔═╡ 4f527a3c-5e0d-4b7d-916f-ff580eb8749e
q33A = matrixdepot("pascal", 4)

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
]

# ╔═╡ 9254e969-0f77-4cdb-9198-618fee92d1db
md"""
## Q35
"""

# ╔═╡ Cell order:
# ╠═8993d71d-9c83-4b2b-87be-fe5afc2fd359
# ╠═45411a2e-449a-4475-838b-187adabfc6b1
# ╠═60bb7f84-877a-4e32-b50e-8c63a5c2b24b
# ╠═111fa6a0-3d84-420e-b546-59075a03a5c2
# ╠═5b4a23b6-17ee-4ce9-814f-cb8dd2489a58
# ╠═49f2380f-af3b-4b20-afde-1b793499cdec
# ╠═2b5d4b30-aeab-11eb-0b6d-5521705ce8f5
# ╠═b18b009f-f2ea-4093-a91c-1bf56c661a86
# ╠═33714c08-90c6-4044-ad80-c5649024f47f
# ╠═a09e4c71-a672-41f8-8989-14f583da0006
# ╠═de645e81-575a-416a-ade9-728ac5f915bd
# ╠═60ada824-3740-40c0-bc55-0946283db92e
# ╠═660eea34-7f99-4214-bccd-1ee0a299405a
# ╠═61948f65-5a96-415d-88f9-d05169b7d775
# ╠═478ce283-3b5c-4a19-8aac-a8e7678928b0
# ╠═c4b032d4-9ab4-4ebd-baed-8fe2ca0a9a64
# ╠═a041250a-6a9a-47d8-9dbe-ee562722a53f
# ╠═daf13d87-4ff6-4cdb-8422-c43eb426ac57
# ╠═263e7145-0eb5-4e43-a5d5-5a714791ef1d
# ╠═3d512797-1045-4c78-9a17-c6aa0d7a80a9
# ╠═672ce5a4-0f6a-48ce-bb54-dad1371ea60c
# ╠═03b9475f-efdf-4897-8797-e8498f343062
# ╠═d2f3b244-dc7f-4a51-99bb-cfeffb64d21b
# ╠═9361ae3c-ab6a-4c9b-8eee-bac81c78eb98
# ╠═c3bfaea3-db75-456a-ae17-9bdfc5cacdef
# ╠═6a133cee-048e-4c26-9c11-be71a22edddc
# ╠═fbe30313-f9e9-4151-b426-df61b56dc28f
# ╠═433b3e27-fd82-4b19-ad90-102618fcc381
# ╠═eddef1e3-ffff-400e-af18-88cfe67a81b2
# ╠═339633aa-1eb7-42b4-b922-6a57fcc57455
# ╠═68a2d75c-04da-4a24-b823-d1ec781289c1
# ╠═269120fd-b9b2-4e90-9c67-c83a0dc7d779
# ╠═9a7b4489-1f5a-4589-9d93-d26fcc470dca
# ╠═23dd8a67-dbf5-4870-b278-d34336dfac30
# ╠═8f77d447-23c2-4832-95f1-c66919e1f3d7
# ╠═208af161-d36c-43b7-adc8-16ad987b8de7
# ╠═58bd585b-ac46-40d2-b51f-10c0f56225d2
# ╠═e4570afd-c1b1-49e5-a4a8-29a092a45c61
# ╠═b4dde61d-44fc-42ed-9e33-4efcbfcc9bc3
# ╠═221f3f2d-442e-48bd-aa2c-5fba5e069844
# ╠═ac3a60fc-4ac4-40a7-9c68-a0bf860dea01
# ╠═44fff8ac-9940-4cb4-9ab2-ee5ea9ced918
# ╠═2fa48fd0-9047-4308-b31d-f0ea5dc4b91a
# ╠═2e6cdf1e-295b-4554-8c13-9737eabbb59a
# ╠═f9710878-110e-4ef7-a3a1-ed47f22d7dcc
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
# ╠═07105032-916f-4caa-850f-fdacd8a1b971
# ╠═c3428cf4-cf48-4b9d-ba3c-da8187018bdd
# ╠═e2c7912a-d980-466a-8843-ae13893839b6
# ╠═6d0194a4-ff29-49dd-8146-1de35f77998b
# ╠═fbe51c6c-5d45-4204-b7e3-d43bb09ecc8b
# ╠═6ae1dd18-b181-4d9d-8246-167cf9535df0
# ╠═7f261229-8220-4be5-980e-0141c3cc31e1
# ╠═86fd8e42-2abc-4b67-9950-1480f7ba636d
# ╠═e7837ea4-f805-45e1-9534-347d253a2de3
# ╠═50aee535-4ec6-467f-ac00-9902b3ead43c
# ╠═99e8e7a2-8591-4f0c-8830-07cc0883b0b6
# ╠═b9172bca-f008-4579-b4c4-d1b7fc1e9a41
# ╠═3235d955-ad52-4afd-8eb4-e7f6eb9dd03f
# ╠═84449c10-438a-44b8-ac1e-dd617296dc8a
# ╠═e4eb0ee4-76eb-4b63-ae81-6efaf90c5338
# ╠═36dd5725-1589-4f29-ae6c-9f4ebe223e83
# ╠═a886d29e-bce2-4316-860d-3063df707f36
# ╠═4447861d-56f6-44ae-91d6-87a393be0a43
# ╠═eb0c5ef2-28d5-431f-8d22-7db0564c24a3
# ╠═6a6ec27c-4cd1-429b-b881-f92552550118
# ╠═63e2463d-fa27-44cf-a2cc-fc53caa019f6
# ╠═13b14998-2bf9-4e5d-b902-d480a5d9711c
# ╠═24ddccd5-8787-4c58-8c62-4dca6787762f
# ╠═03399151-21e9-4e5f-aa6e-a7813ca47a15
# ╠═cb514a8a-123c-4839-a17d-8d1fb41430ab
# ╠═59d66bf5-20d9-4eba-b245-28341aec58c0
# ╠═0be7d6ab-313c-4251-a3fc-014d0ff08d8e
# ╠═8baed9ad-8824-41f8-ab8a-c2acf60edb46
# ╠═0a8fad3e-759e-4b8b-a3a0-466574feb9c9
# ╠═4f527a3c-5e0d-4b7d-916f-ff580eb8749e
# ╠═f8653bd3-d2f0-43fd-9d64-381aa7a5c5aa
# ╠═dfe40c17-4f5d-4a63-b6d1-cb4d6f4d178d
# ╠═9254e969-0f77-4cdb-9198-618fee92d1db
