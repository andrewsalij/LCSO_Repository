import numpy as np

e_xy = np.array([[1j,0],[0,2j]])

vals,vecs= np.linalg.eig(e_xy)

vals2,vecs2 = np.linalg.eig(np.real(e_xy))
vals3,vecs3 = np.linalg.eig(np.imag(e_xy))