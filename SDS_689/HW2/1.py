import numpy as np   

# generate the dataset
N_data = 200 # number of points for data
N_test = 689 # number of points for testing
k_tot = 4 # number of elements in sinusoidal sum
X = np.linspace(-np.pi,np.pi,N_data) # create grid ofN_data points in x

# define function to evaluate f on grid points in X
def func(xq):
    yq = np.zeros(len(xq))
    # separate X into segments below and above zero
    x_below = xq[xq < 0]
    x_above = xq[xq >= 0]
    yq[xq < 0] = 5
    for k in (np.arange(k_tot) + 1):
        yq[xq < 0] += np.sin(k * x_below)
    yq[xq >= 0] = np.cos(10 * x_above)
    return yq
y = func(X)

# define L2 relative error
def L2_rel_error(approx, actual):
    top = np.sum(np.square(approx - actual))
    bottom = np.sum(np.square(actual))
    return np.sqrt(top / bottom)

# plot function
import matplotlib.pyplot as plt
plt.plot(X, y)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of discontinuous oscillatory function from Problem 1')
plt.savefig('1_plot.pdf')
plt.close()

# train neural network
import torch
import torch.nn as nn
import torch.optim as optim

if torch.cuda.is_available():
  print("CUDA available. Using GPU acceleration.")
  device = "cuda"
else:
  print("CUDA is NOT available. Using CPU for training.")
  device = "cpu"

# convert input and output data to PyTorch tensors
X = torch.tensor(X, dtype=torch.float32, device=device).reshape(-1, 1)
y = torch.tensor(y, dtype=torch.float32, device=device).reshape(-1, 1)

# testing domain
X_test = torch.tensor(np.linspace(-np.pi,np.pi,N_test), dtype=torch.float32, 
                      device=device).reshape(-1, 1)

widths = [10, 30, 100, 300, 1000] # list of widths to test
loss_fn = nn.MSELoss()  # use MSE loss
n_epochs = 5e4 # number of epochs
lr = 0.001 # learning rate for Adam optimizer 
for w in widths:
    # define shallow ReLU network architecture with width w 
    model = nn.Sequential(
        nn.Linear(1, w),
        nn.ReLU(),
        nn.Linear(w, 1)).to(device)
    # use Adam optimizer with learning rate lr
    optimizer = optim.Adam(model.parameters(), lr=lr) 
    for epoch in (np.arange(n_epochs) + 1):
        y_pred = model(X) # forward pass
        loss = loss_fn(y_pred, y) # compute loss
        optimizer.zero_grad() # zero out gradient
        loss.backward() # backpropagation
        optimizer.step() # adjust parameters
    print(f'Finished training with width {w}, latest loss {loss}')
    # testing
    y_pred = torch.flatten(model(X_test)).detach().cpu().numpy()
    X_testq = torch.flatten(X_test).cpu().numpy()
    y_actual = func(X_testq)
    plt.plot(X_testq, y_pred, label='Prediction')
    plt.plot(X_testq, y_actual, label='Actual')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Prediction vs. actual from shallow ReLU network \n '+
            'Network width: '+str(w)+', Number of epochs: '+str(n_epochs))
    plt.savefig('1_testing_'+str(w)+'.pdf')
    plt.close()
    # return L2 relative error for width w
    print(f'L2 relative error for width {w}: {L2_rel_error(y_pred, y_actual)}') 
        

