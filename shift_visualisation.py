import numpy as np 
import matplotlib.pyplot as plt

shifts0 = np.load('tmp_1005/flow16.npy')
shifts1 = np.load('tmp_1005/flow20.npy')
shifts2 = np.load('tmp_1005/flow24.npy')
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