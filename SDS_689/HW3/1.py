# download digits datasets

# the CSVs are in the following format:
# -first column is the digit in question
# -rest of the columns contain the image data
# -rows are samples

import pandas as pd
import numpy as np
import torch

if torch.cuda.is_available():
  print("CUDA available. Using GPU acceleration.")
  device = "cuda"
else:
  print("CUDA is NOT available. Using CPU for training.")
  device = "cpu"

# load digit dataset for TRAINING
# load directly as array of ints, skip header
train_data = pd.read_csv('mnist_train.csv', skiprows=1, header=None, dtype=int).values
# load labels as a list of unit vectors 
# (0 if wrong label, 1 if correct label, for each sample)
train_labels = torch.tensor(train_data[:,0], device=device)
train_labels_tensor = torch.tensor(np.arange(10) == train_data[:,0,None], 
                                   dtype=torch.float32, device=device)
n_train = train_data[:,0].size # number of samples in training
N = int(np.sqrt(train_data[0].size - 1)) # length of 2D image
# reshape each 1D array of sample data of size 784 to a 2D image of size (28, 28, 1)
train_data = np.expand_dims(train_data[:,1:].reshape(n_train, N, N), axis=1)
# normalize train_data, convert to PyTorch tensor
train_data = torch.tensor(train_data / np.max(train_data), 
                          dtype=torch.float32, device=device) 
        
# load digit dataset for TESTING
# load directly as array of ints, skip header
test_data = pd.read_csv('mnist_test.csv', skiprows=1, header=None, dtype=int).values
n_test = test_data[:,0].size
test_labels = torch.tensor(test_data[:,0], device=device)
test_data = np.expand_dims(test_data[:,1:].reshape(n_test, N, N), axis=1)
# normalize train_data, convert to PyTorch tensor
test_data = torch.tensor(test_data / np.max(test_data), 
                         dtype=torch.float32, device=device)

# train neural network
import torch.nn as nn
import torch.optim as optim

# hyperparameters
loss_fn = nn.CrossEntropyLoss()  # use cross entropy loss
n_epochs = 50 # number of epochs
lr = 0.001 # learning rate for Adam optimizer 
batches = 1e2 # number of batches

# use Sequential PyTorch class to construct LeNet-5 architecture
# with some differences:
# -replaced Gaussian filter with ReLU
# -added batch normalization between convolution and ReLU layers

class CNN(nn.Module):
    def __init__(self, channels, classes):
        super(CNN, self).__init__() # initialize parent constructor
        
	    # conv => batchnorm => ReLU => pooling (first)
        self.layer1 = nn.Sequential(
            nn.Conv2d(channels, 6, kernel_size=(5, 5), padding=2),
            nn.BatchNorm2d(6),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(2, 2), stride=2))
        
		# conv => batchnorm => ReLU => pooling (second)
        self.layer2 = nn.Sequential(
            nn.Conv2d(6, 16, kernel_size=(5, 5), padding=0),
            nn.BatchNorm2d(16),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(2, 2), stride=2))
        
        # series of FC NNs
        self.layer3 = nn.Sequential(
            nn.Linear(16 * 5 * 5, 120),
            nn.ReLU(),
            nn.Linear(120, 84),
            nn.ReLU(),
            nn.Linear(84, classes)
        )
    
    def forward(self, x):
        output = self.layer1(x)
        output = self.layer2(output)
        output = torch.flatten(output,start_dim=1)
        output = self.layer3(output)
        return output

model = CNN(channels=1, classes=10).to(device) # define network
# use Adam optimizer with learning rate lr
optimizer = optim.Adam(model.parameters(), lr=lr) 
batch_size = int(np.floor(n_train / batches)) # batch size
from random import shuffle 
for epoch in (np.arange(n_epochs) + 1):
    inds = list(range(n_train))
    shuffle(inds)
    running_loss = 0
    for i in np.arange(batches): # loop through batches
        # randomly sample indices for each batch without replacement 
        batch_ind = inds[-batch_size:] 
        del inds[-batch_size:]
        pred = model(train_data[[batch_ind]]) # forward pass
        loss = loss_fn(pred, train_labels_tensor[[batch_ind]]) # compute loss
        running_loss += loss * len(batch_ind)
        optimizer.zero_grad() # zero out gradient
        loss.backward() # backpropagation
        optimizer.step() # adjust parameters
    print(f'Cross-entropy loss for epoch {epoch}: {running_loss / n_train}')
    predicted_labels = torch.argmax(model(train_data), dim=1)
    frac_incorrect = torch.count_nonzero(predicted_labels - train_labels) / n_train
    print(f'Training error for epoch {epoch}: {frac_incorrect}')
    
# test network on test dataset, get predicted labels
predicted_labels = torch.argmax(model(test_data), dim=1)
# test fraction of predicted labels correct
frac_incorrect = torch.count_nonzero(predicted_labels - test_labels) / n_test
print(f'Testing error: {frac_incorrect}')