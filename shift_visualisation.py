import numpy as np 
import matplotlib.pyplot as plt

shifts0 = np.load('1127_1005n/tmp_1005_/flow32.npy')
shifts1 = np.load('1127_1005n/tmp_1005_/flow48.npy')
shifts2 = np.load('1127_1005n/tmp_1005_/flow64.npy')
# shifts0 = np.load('1150_0/tmp_1005_1e-06/flow8.npy')
# shifts1 = np.load('1150_0/tmp_1005_1e-06/flow12.npy')
# shifts2 = np.load('1150_0/tmp_1005_1e-06/flow16.npy')

plt.subplot(2,1,1)
plt.title('vertical shifts')
plt.plot(shifts0[:,0],'r.')
plt.plot(shifts1[:,0],'g.')
plt.plot(shifts2[:,0],'b.')
plt.subplot(2,1,2)
plt.title('horizontal shifts')
plt.plot(shifts0[:,1],'r.')
plt.plot(shifts1[:,1],'g.')
plt.plot(shifts2[:,1],'b.')
plt.show()