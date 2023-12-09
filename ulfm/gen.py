import numpy as np

def rand_inv(n):
  m = np.random.rand(n, n)
  mx = np.sum(np.abs(m), axis=1)
  np.fill_diagonal(m, mx + 1)
  return m

def with_rank(n, m ,r):
  a = rand_inv(n)
  b = rand_inv(m)

  rr = np.zeros((n, m))
  for i in range(r):
    rr[i,i] = 1

  res = ((a.dot(rr)).dot(b))
  print("rank is", np.linalg.matrix_rank(res))
  return res.flatten()

N = 3*100
M = 200

line_len = 47
rank = 123

arr = with_rank(N, M, rank)

with open("rank.mt", "w") as f:
  for e in arr:
    f.write(f'{e}\n')

with open("rank.inc", "w") as f:
  f.write(
      f"""    parameter (N = {N}, M = {M})\n    real, dimension(N, M):: A\n"""
  )

print("rank is" , rank)