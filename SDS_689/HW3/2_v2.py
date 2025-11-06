a = 0 # start of training domain
b = 10 # end of training domain/start of testing domain
c = 12 # end of testing domain 
# (training is on [a,b], testing is on [b,c])
dx = 0.05 # interval size
mu = 0 # mean of Gaussian noise
sigma = 0.02 # standard deviation of Gaussian noise
N = 10 # lookback period; number of prior consecutive samples used to predict the next
n_epochs = 100 # number of epochs for training
w = 50 # width of hidden dense layer
lr = 0.001 # learning rate

# generate training dataset
import numpy as np
from keras import Sequential
from keras.layers import LSTM, Dense
import matplotlib.pyplot as plt

## Data Preperation
def seq_data(x, seq_length):
    X = []
    Y = []
    l = len(x)
    
    for i in range(l):
        end_id = i + seq_length
        if end_id > len(x) - 1:
            print("end id could not be bigger than series length")
            break    
        X.append(x[i:end_id])
        Y.append(x[end_id])
        
    return np.array(X), np.array(Y)
    
x_train = np.linspace(a, b, int((b-a) / dx))
y_train = np.cos(2*np.pi*x_train)
# add noise to training data
y_train += np.random.normal(loc=mu, scale=sigma, size=np.shape(y_train)) 
x_test = np.linspace(b, c, int((c-b) / dx))
y_test = np.cos(2*np.pi*x_test)

x_train_plot = np.copy(x_train)
y_train_plot = np.copy(y_train)
x_test_plot = np.copy(x_test)
y_test_plot = np.copy(y_test)

x_train, y_train = seq_data(y_train, N) 
x_test, y_test = seq_data(y_test, N) 

train_shape = x_train.shape
test_shape = x_test.shape

x_train = x_train.reshape((train_shape[0], train_shape[1], 1))
x_test = x_test.reshape((test_shape[0], test_shape[1], 1))

from keras.optimizers import Adam
optimizer = Adam(lr=lr)
model = Sequential()
model.add(LSTM(w, input_shape = (N, 1)))
model.add(Dense(1))
model.compile(optimizer=optimizer, loss='mean_squared_error')

h = model.fit(x_train, y_train, epochs=n_epochs, verbose=1)
# add N rows to predict N steps after the end of training dataset
to_add = np.zeros((N, N, 1))
for i in np.arange(N):
    if i == 0:
        to_add[i] = y_train[-N:, None]
    else:
        to_add[i] = y_train[(-N-i):(-i), None]
to_add[:-1] = to_add[1:]
to_add[-1] = np.insert(to_add[-2,:-1], 0, np.cos(2*np.pi*b))[:,None]
# test model to predict function on [10,12]
y_predict = model.predict(np.concatenate((to_add, x_test)))[:,0] 

# define L2 relative error
def L2_rel_error(approx, actual):
    top = np.sum(np.square(approx - actual))
    bottom = np.sum(np.square(actual))
    return np.sqrt(top / bottom)

# return L2 relative testing error
print(f'L2 relative testing error: {L2_rel_error(y_predict, y_test_plot)}') 

plt.plot(x_train_plot, y_train_plot, label='Training data')
plt.plot(x_test_plot, y_test_plot, label='Testing data')
plt.plot(x_test_plot, y_predict, label='Prediction via LSTM')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.show()