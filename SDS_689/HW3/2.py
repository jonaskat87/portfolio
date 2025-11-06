a = 0 # start of training domain
b = 10 # end of training domain/start of testing domain
c = 12 # end of testing domain 
# (training is on [a,b], testing is on [b,c])
dx = 0.05 # interval size
mu = 0 # mean of Gaussian noise
sigma = 0.02 # standard deviation of Gaussian noise
n_epochs = 100 # number of epochs for training
w = int(input()) # width of hidden dense layer (taken as command line input)
lr = 0.001 # learning rate 

# generate training dataset
import numpy as np
from keras import Sequential
from keras.layers import LSTM, Dense
import matplotlib.pyplot as plt

x_train = np.linspace(a, b, int((b-a) / dx))
y_train = np.cos(2*np.pi*x_train)
# add noise to training data
y_train += np.random.normal(loc=mu, scale=sigma, size=np.shape(y_train)) 
x_test = np.linspace(b, c, int((c-b) / dx))
y_test = np.cos(2*np.pi*x_test)

input_train = y_train.reshape((len(y_train),1,1))
input_test = y_test.reshape((len(y_test),1,1))

from keras.optimizers import Adam
optimizer = Adam(learning_rate=lr)
model = Sequential()
model.add(LSTM(w, input_shape = (1, 1)))
model.add(Dense(1))
model.compile(optimizer=optimizer, loss='MSE')

h = model.fit(input_train, y_train, epochs=n_epochs, verbose=1)
y_pred = model.predict(input_test).reshape(1,-1)

# define L2 relative error
def L2_rel_error(approx, actual):
    top = np.sum(np.square(approx - actual))
    bottom = np.sum(np.square(actual))
    return np.sqrt(top / bottom)

# return L2 relative testing error
print(f'L2 relative testing error: {L2_rel_error(y_pred, y_test)}') 

plt.title(f'LSTM prediction of $\cos(2\pi x)$ with width={w}')
plt.plot(x_train, y_train, label='Training data')
plt.plot(x_test, y_test, label='Testing data')
plt.plot(x_test, y_pred[0], label='Prediction via LSTM')
plt.xlabel('$x$')
plt.ylabel('$f(x)=\cos(2\pi x)$')
plt.legend()
plt.savefig(f'2_{w}.pdf')
plt.close()