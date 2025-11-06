# generate training dataset
import numpy as np
a = 0 # start of training domain
b = 10 # end of training domain/start of testing domain
c = 12 # end of testing domain 
# (training is on [a,b], testing is on [b,c])
dx = 0.05 # interval size
mu = 0 # mean of Gaussian noise
sigma = 0.02 # standard deviation of Gaussian noise
N = 10 # lookback period; number of prior consecutive samples used to predict the next

xq_training = np.linspace(a, b, int((b-a) / dx))
xq_testing = np.linspace(b, c, int((c-b) / dx))
fq_training = np.cos(2*np.pi*xq_training) # evaluate f(x) on grid for training
fq_training += np.random.normal(loc=mu, scale=sigma, size=np.shape(fq_training)) # add noise to training data
fq_testing = np.cos(2*np.pi*xq_testing) # evaluate f(x) on grid (ground truth) for testing

# split fq_training into inputs and outputs, assuming LSTM uses a lookback period of N
# rows are samples, each column is a dimension. This can be stored using a structure called a Hankel matrix.

import torch

if torch.cuda.is_available():
  print("CUDA available. Using GPU acceleration.")
  device = "cuda"
else:
  print("CUDA is NOT available. Using CPU for training.")
  device = "cpu"

from scipy.linalg import hankel
hankel_training = hankel(fq_training)
input_training = torch.tensor(hankel_training[:-N,:N], 
                              dtype=torch.float32, device=device)
output_training = torch.tensor(hankel_training[:-N,N].reshape(-1, 1), 
                               dtype=torch.float32, device=device)
n = input_training.size(dim=0)

# define L2 relative error
def L2_rel_error(approx, actual):
    top = np.sum(np.square(approx - actual))
    bottom = np.sum(np.square(actual))
    return np.sqrt(top / bottom)

# train RNN (LSTM)
import torch.nn as nn
import torch.optim as optim

# define hyperparameters
loss_fn = nn.MSELoss()  # use MSE loss
n_epochs = 20000 # number of epochs
lr = 0.01 # learning rate for Adam optimizer
h = 20 # size of hidden layer
k = 100 # print every kth epoch
batches = 5 # number of batches

# define model architecture
class recurrent(nn.Module):
    def __init__(self, lookback_period, hidden_dim):
        super(recurrent, self).__init__() # initialize parent constructor
        self.lstm1 = nn.LSTM(lookback_period, hidden_dim, batch_first=True) 
        self.lstm2 = nn.LSTM(hidden_dim, hidden_dim, batch_first=True) 
        self.hidden2 = nn.Linear(hidden_dim, 1)
    
    def forward(self, x):
        output, _ = self.lstm1(x)
        output, _ = self.lstm2(output[:, -1, :])
        output = self.hidden2(output)
        return output
    
model = recurrent(lookback_period=N, hidden_dim=h).to(device) # initialize NN
optimizer = optim.Adam(model.parameters(), lr=lr) # use Adam optimizer with learning rate lr
batch_size = int(np.floor(n / batches)) # batch size
from random import shuffle 
for epoch in (np.arange(n_epochs) + 1):
    inds = list(range(n))
    shuffle(inds)
    running_loss = 0
    for i in np.arange(batches): # loop through batches
        # randomly sample indices for each batch without replacement 
        batch_ind = inds[-batch_size:] 
        del inds[-batch_size:]
        y_pred = model(input_training[[batch_ind]])
        loss = loss_fn(y_pred, output_training[[batch_ind]])
        running_loss += loss * len(batch_ind)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    if epoch%k == 0: # print every kth epoch
        print(f'MSE loss for epoch {epoch}: {running_loss / n}')
        y_pred_training = model(input_training).detach().cpu().numpy()
        err = L2_rel_error(y_pred_training, output_training.cpu().numpy())
        print(f'L2 relative training error: {err}') 

# compute prediction on [b,c] recursively 
# initialize prediction
y_pred = torch.empty(np.shape(fq_testing), dtype=torch.float32, device=device) 
# initialize input vector
input_testing = torch.tensor(fq_training[None,-N:], dtype=torch.float32, device=device)
y_pred[0] = torch.tensor(fq_training[None,-1], dtype=torch.float32, device=device)
for i in np.arange(fq_testing.size - 1):
    i += 1
    y_pred[i] = model(input_testing) # predict next timestep using LSTM
    input_testing = torch.cat((input_testing[:,1:],y_pred[i].reshape(1,1)),1)
    
# return L2 relative testing error
y_pred_testing = y_pred.detach().cpu().numpy()
print(f'L2 relative testing error: {L2_rel_error(y_pred_testing, fq_testing)}') 

import matplotlib.pyplot as plt
plt.plot(xq_training, fq_training, label='Training data')
plt.plot(xq_testing, fq_testing, label='Testing data')
xq = np.concatenate((xq_training[N:-1], xq_testing))
y_pred = np.concatenate((y_pred_training[:-1,0], y_pred_testing))
plt.plot(xq, y_pred, label='Prediction via LSTM')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.show()