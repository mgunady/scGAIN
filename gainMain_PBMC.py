'''
Written by Mohamed Gunady
Date: Dec 2018 (Revised Dec 2nd 2019)
scGAIN: Single Cell RNA-seq Data Imputation using Generative Adversarial Networks
Paper Link: https://www.biorxiv.org/content/10.1101/837302v1
Contact: mgunady@cs.umd.edu
-------------
Based on code written by Jinsung Yoon
Date: Jul 9th 2018 (Revised Oct 19th 2018)
Generative Adversarial Imputation Networks (GAIN) Implementation on MNIST
Reference: J. Yoon, J. Jordon, M. van der Schaar, "GAIN: Missing Data Imputation using Generative Adversarial Nets," ICML, 2018.
Paper Link: http://medianetlab.ee.ucla.edu/papers/ICML_GAIN.pdf
Appendix Link: http://medianetlab.ee.ucla.edu/papers/ICML_GAIN_Supp.pdf
Contact: jsyoon0823@g.ucla.edu
'''

#%% Packages
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data
import numpy as np
import os
from tqdm import tqdm

import tflib as lib
import tflib.plot
import math

import pandas as pd

def calcNMSE(truth, pred):
	mse = np.sum(np.abs(truth-pred), axis=1)
	truthNorm = np.sum(truth, axis=1)
	print(mse.shape, truthNorm.shape, min(mse), min(truthNorm))
	nmse = mse/truthNorm
	print(np.mean(nmse), np.median(nmse))
	return(np.mean(nmse), np.median(nmse))

def preprocessData(X, maxZerosInCell, maxZerosInGene):
	cellSums = (X==0).sum(axis=1)
	selectedCells = cellSums <= maxZerosInCell
	geneSums = (X==0).sum(axis=0)
	selectedGenes = geneSums <= maxZerosInGene
	#print(cellSums, geneSums, maxZerosInCell, maxZerosInGene)
	#print(geneSums, np.sum(selectedGenes))
	selectedCellsIdxs = np.array([i for i, v in enumerate(selectedCells) if v])
	selectedGenesIdxs = np.array([i for i, v in enumerate(selectedGenes) if v])
	#print(selectedCellsIdxs, selectedGenesIdxs)
	X_f = X[selectedCellsIdxs[:, None], selectedGenesIdxs]
	X_f_log = np.log(X_f+1)
	#print("==============")
	#print(X[:5, :5])
	#print(X_f[:5, :5])
	#print(X_f_log[:5, :5])
	maxPerCell = X_f_log.max(axis=1)
	print(np.min(X), np.max(X), np.min(maxPerCell), np.max(maxPerCell))
	#print(maxPerCell[:5], len(maxPerCell))
	X_f_log_norm = X_f_log/maxPerCell[:,None]
	return(X_f_log_norm, selectedGenes, selectedCells, maxPerCell)

def preprocessData2(X, maxZerosInCell, maxZerosInGene):
	cellSums = (X==0).sum(axis=1)
	selectedCells = cellSums <= maxZerosInCell
	geneSums = (X==0).sum(axis=0)
	selectedGenes = geneSums <= maxZerosInGene
	#print(cellSums, geneSums, maxZerosInCell, maxZerosInGene)
	#print(geneSums, np.sum(selectedGenes))
	selectedCellsIdxs = np.array([i for i, v in enumerate(selectedCells) if v])
	selectedGenesIdxs = np.array([i for i, v in enumerate(selectedGenes) if v])
	#print(selectedCellsIdxs, selectedGenesIdxs)
	X_f = X[selectedCellsIdxs[:, None], selectedGenesIdxs]
	
	s = X_f.sum(axis=1)
	s_m = np.median(s)
	s_ = s/s_m
	X_f_norm = X_f/s_[:,None]
	X_f_norm_log = np.log(X_f_norm+1)
	#print("==============")
	#print(X[:5, :5])
	#print(X_f[:5, :5])
	#print(X_f_log[:5, :5])
	return(X_f_norm_log, selectedGenes, selectedCells, s_)


def nonZeroMean(v):
	meu = 0 if np.count_nonzero(v) == 0 else np.median(v[v!=0]).round(decimals=2)
	return(meu)

def getMask(X, alpha, sparsity=101):
	geneAvgs = pd.DataFrame(X).apply(nonZeroMean)
	geneSums = (X==0).sum(axis=0)
	maxNZerosInGene = sparsity*X.shape[0]/100
	#print(geneAvgs)
	#print(geneAvgs)
	mask = np.ones(X.shape)
	#X_bar = mask*geneAvgs[None,:]
	mask[(geneAvgs[None,:] > alpha) & (geneSums <= maxNZerosInGene) & (X == 0)] = 0
	return(mask, geneAvgs)

def transformBack(X_, X, M):
	X_ = M*X+(1-M)*X_
	X_t = np.transpose(X_)
	return(X_t)

def transformBackAll(X_, X, M, filteredGenes, filteredCells, maxPerCell):
	X_ = M*X+(1-M)*X_
	print(np.min(X_), np.max(X_))
	X_ = X_*maxPerCell[:,None]
	print(np.min(X_), np.max(X_))
	X_ = np.exp(X_)-1
	print(np.min(X_), np.max(X_))
	X_t = np.transpose(X_)
	return(X_t)


# Rows are genes, Cols are cells
data_suffix = 'PBMC'#'dropout_index_5_seed_20000'
out_suffix = 'PBMC'#'5_20000'
data = pd.read_csv('simulation_data/simulation_data_'+data_suffix+'.csv', delimiter=',', header=None)
data_full = pd.read_csv('simulation_data/simulation_data_'+data_suffix+'_logn_true.csv', delimiter=',', header=None)
data_full = data_full.T.to_numpy()
#data = data.loc[1:10,1:6].T
#print(data.to_numpy)
data = data.to_numpy()#[:1000,:]
print("Data with ", data.shape[0], " cells, and ", data.shape[1], " genes")
#print(data[:5, :5])

maxZerosInCell = 95*data.shape[1]/100
maxZerosInGene = 95*data.shape[0]/100
small_zero_thre = 2#0.6#0.7
data_f, filteredGenes, filteredCells, maxPerCell = preprocessData2(data, maxZerosInCell, maxZerosInGene)


selectedCellsIdxs = np.array([i for i, v in enumerate(filteredCells) if v])
selectedGenesIdxs = np.array([i for i, v in enumerate(filteredGenes) if v])
print(len(selectedCellsIdxs), len(selectedGenesIdxs), selectedCellsIdxs, selectedGenesIdxs)
data_full = data_full[selectedGenesIdxs[:, None], selectedCellsIdxs]
np.savetxt('imputation_gain_data/'+out_suffix+'_selectedGenes.csv', selectedGenesIdxs+1 , delimiter=',', fmt="%d")



#print(data_f, filteredGenes, filteredCells, maxPerCell)
mask, geneAvgs = getMask(data_f, small_zero_thre, 101)
#data_f = data_f.to_numpy()

print("Impute Matrix After Preprocessing ", data_f.shape)
#print(mask.shape)
print(data[:5, :5])
print(data_f[:5, :5])
print(mask[:5, :5])
np.savetxt('imputation_gain_data/'+out_suffix+'_logn.csv', data_f, delimiter=',', fmt="%f")
#np.savetxt('imputation_gain_data/'+out_suffix+'_mask.csv', mask, delimiter=',', fmt="%f")
#mask2 = pd.read_csv('simulation_data/mask.csv', delimiter=',', header=None).T.to_numpy()
print("Mask Dim ", mask.shape)
#print(mask[:5, :5])
#res = mask2*2+mask
#np.savetxt('imputation_gain_data/'+out_suffix+'_mask_diff.csv', res, delimiter=',', fmt="%i")
#np.savetxt('imputation_gain_data/'+out_suffix+'_mask_avgs.csv', geneAvgs, delimiter=',', fmt="%.3f")
#idxs = np.where(res[1,]==2)
#print(geneAvgs[idxs[0]], data_f[1, idxs[0]])
#exit(0)

#%% System Parameters
# 1. Mini batch size
mb_size = 128
# 3. Hint rate
p_hint = 0.9
# 4. Loss Hyperparameters
alpha = 5
# 5. Imput Dim (Fixed)
Dim = data_f.shape[1]

# Mask Vector and Hint Vector Generation
def sample_M(m, n, p):
	A = np.random.uniform(0., 1., size = [m, n])
	B = A > p
	C = 1.*B
	#C[:,7150:7155] = 0
	return C

def sample_M_bias(m, n, p, probs):
    #probs = probs/sum(probs)
    num = int(p*n)
    l = np.array([np.random.choice(n, num, False, probs) for i in range(m)])
    rows = np.repeat(range(m), num)
    cols = l.reshape(-1)
    #print(l, rows, cols)

    mask = np.ones((m, n))
    mask[rows, cols] = 0
    return(mask)


######################
## TensorFlow
######################

#%% Necessary Functions
# 1. Xavier Initialization Definition
def xavier_init(size):
    in_dim = size[0]
    xavier_stddev = 1. / tf.sqrt(in_dim / 2.)
    return tf.random_normal(shape = size, stddev = xavier_stddev)
    
# 2. Plot (4 x 4 subfigures)
def plot(samples):
    fig = plt.figure(figsize = (5,5))
    gs = gridspec.GridSpec(5,5)
    gs.update(wspace=0.05, hspace=0.05)
    
    for i, sample in enumerate(samples):
        ax = plt.subplot(gs[i])
        plt.axis('off')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
        plt.imshow(sample.reshape(298,11), cmap='Greys_r')
        
    return fig
   
'''
GAIN Consists of 3 Components
- Generator
- Discriminator
- Hint Mechanism
'''   
   
#%% GAIN Architecture   
   
#%% 1. Input Placeholders
# 1.1. Data Vector
X = tf.placeholder(tf.float32, shape = [None, Dim])
# 1.2. Mask Vector 
M = tf.placeholder(tf.float32, shape = [None, Dim])
# 1.3. Hint vector
H = tf.placeholder(tf.float32, shape = [None, Dim])
# 1.4. Random Noise Vector
Z = tf.placeholder(tf.float32, shape = [None, Dim])
NZ = tf.placeholder(tf.float32, shape = [None, Dim])

#%% 2. Discriminator
D_W1 = tf.Variable(xavier_init([Dim*2, 256]))     # Data + Hint as inputs
D_b1 = tf.Variable(tf.zeros(shape = [256]))

D_W2 = tf.Variable(xavier_init([256, 128]))
D_b2 = tf.Variable(tf.zeros(shape = [128]))

D_W3 = tf.Variable(xavier_init([128, Dim]))
D_b3 = tf.Variable(tf.zeros(shape = [Dim]))       # Output is multi-variate

theta_D = [D_W1, D_W2, D_W3, D_b1, D_b2, D_b3]

#%% 3. Generator
G_W1 = tf.Variable(xavier_init([Dim*2, 256]))     # Data + Mask as inputs (Random Noises are in Missing Components)
G_b1 = tf.Variable(tf.zeros(shape = [256]))

EMB_SIZE=128
G_W2 = tf.Variable(xavier_init([256, EMB_SIZE]))
G_b2 = tf.Variable(tf.zeros(shape = [EMB_SIZE]))

G_W3 = tf.Variable(xavier_init([EMB_SIZE, Dim]))
G_b3 = tf.Variable(tf.zeros(shape = [Dim]))

theta_G = [G_W1, G_W2, G_W3, G_b1, G_b2, G_b3]

#%% GAIN Function

#%% 1. Generator
def generator(x,z,m):
    inp = m * x + (1-m) * z  # Fill in random noise on the missing values
    inputs = tf.concat(axis = 1, values = [inp,m])  # Mask + Data Concatenate
    G_h1 = tf.nn.relu(tf.matmul(inputs, G_W1) + G_b1)
    G_h2 = tf.nn.relu(tf.matmul(G_h1, G_W2) + G_b2)
    #G_prob = tf.nn.sigmoid(tf.matmul(G_h2, G_W3) + G_b3) # [0,1] normalized Output
    G_prob = tf.nn.relu(tf.matmul(G_h2, G_W3) + G_b3) # [0,1] normalized Output
    
    return G_prob, G_h2
    
#%% 2. Discriminator
def discriminator(x, m, g, h):
    inp = m * x + (1-m) * g  # Replace missing values to the imputed values
    inputs = tf.concat(axis = 1, values = [inp,h])  # Hint + Data Concatenate
    D_h1 = tf.nn.relu(tf.matmul(inputs, D_W1) + D_b1)
    D_h2 = tf.nn.relu(tf.matmul(D_h1, D_W2) + D_b2)
    D_logit = tf.matmul(D_h2, D_W3) + D_b3
    D_prob = tf.nn.sigmoid(D_logit)  # [0,1] Probability Output
    
    return D_prob

#%% 3. Others
# Random sample generator for Z
def sample_Z(m, n):
    return np.random.uniform(0., 1., size = [m, n])        

def sample_idx(m, n):
    A = np.random.permutation(m)
    idx = A[:n]
    return idx

#%% Structure
G_sample, G_embed = generator(X,Z,M)
D_prob = discriminator(X, M, G_sample, H)

#%% Loss
NZC = 1#NZ#+1-NZ
D_loss1 = -tf.reduce_mean(M * NZC * tf.log(D_prob + 1e-8) + (1-M) * NZC * tf.log(1. - D_prob + 1e-8)) * 2
G_loss1 = -tf.reduce_mean((1-M)  * NZC * tf.log(D_prob + 1e-8)) / tf.reduce_mean(tf.maximum(1-M * NZC, 1) )
MSE_train_loss = tf.reduce_mean( NZC * (M * X - M * G_sample)**2) / tf.reduce_mean(tf.maximum(M * NZC, 1) )

D_loss = D_loss1
G_loss = G_loss1  + alpha * MSE_train_loss 

#%% MSE Performance metric
MSE_test_loss = tf.reduce_mean(((1-M) * NZC  * X - (1-M) * NZC * G_sample)**2) / tf.reduce_mean(tf.maximum((1-M) * NZC, 1) )

#%% Solver
D_solver = tf.train.AdamOptimizer().minimize(D_loss, var_list=theta_D)
G_solver = tf.train.AdamOptimizer().minimize(G_loss, var_list=theta_G)

# Sessions
os.environ["CUDA_VISIBLE_DEVICES"]="0"
saver = tf.train.Saver()
sess = tf.Session()
sess.run(tf.global_variables_initializer())

#%%
# Output Initialization
if not os.path.exists('imputation_gain_data/'):
    os.makedirs('imputation_gain_data/')
    
# Iteration Initialization
trainX = data_f
trainM = mask
testX = trainX[:,:]
testM = trainM[:,:]
Train_No = trainX.shape[0]
geneProbs = np.sum(trainM)

#print([nonZeroMean(trainX[:,i]) for i in range(7150, 7155)])
print("NZeros=", np.count_nonzero(trainX)/(trainX.shape[0]*trainX.shape[1]*1.0), np.count_nonzero(testX)/(trainX.shape[0]*trainX.shape[1]*1.0))
print("NZeros=", np.count_nonzero(1-testM)/(trainX.shape[0]*trainX.shape[1]*1.0), np.count_nonzero(trainX*(1-testM))*1.0/np.count_nonzero(1-testM))

cutoff_i = 0
cutoffs =  [2,    1,     0.75,    0.5, 	  0.5, 0.5]
sparsity = [101,  101,    101,	  60,	 85, 90]
maxIters = [2000, 10000, 25000, 70000,	70000, 70000]
maxIters = [1000, 3000, 6000, 70000,	70000, 70000]
percent_nonzero = 0#np.count_nonzero(trainX*(1-np.transpose(testM)))*1.0/np.count_nonzero(1-testM)

#%% Start Iterations
for it in tqdm(range(70020)):    
    
    #%% Inputs
    mb_idx = sample_idx(Train_No, mb_size)
    X_mb = trainX[mb_idx,:]  
    Z_mb = sample_Z(mb_size, Dim) 
    M_mb = trainM[mb_idx,:]  
    H_mb1 = sample_M(mb_size, Dim, 1-p_hint)
    H_mb = M_mb * H_mb1 + 0.5 * (1-H_mb1)
    
    
    New_X_mb = M_mb * X_mb + (1-M_mb) * Z_mb  # Missing Data Introduce
    
    _, D_loss_curr = sess.run([D_solver, D_loss1], feed_dict = {X: X_mb, M: M_mb, Z: New_X_mb, H: H_mb, NZ:(X_mb > 0)*1})
    _, G_loss_curr, MSE_train_loss_curr, MSE_test_loss_curr = sess.run([G_solver, G_loss1, MSE_train_loss, MSE_test_loss],
                                                                       feed_dict = {X: X_mb, M: M_mb, Z: New_X_mb, H: H_mb, NZ:(X_mb > 0)*1})

    #%% Intermediate Losses
    if it % 100000 == 0:
        print('Iter: {}'.format(it))
        print('Train_loss: {:.4}'.format(MSE_train_loss_curr))
        print('Test_loss: {:.4}'.format(MSE_test_loss_curr))
        print()
        
    #%% Output figure
    if it % 1000 == 0 and it > 1:
    #if it == 3000 or it == 25000 or it == 50000 or it == 70000:
        preds, MSE_train_loss_curr, MSE_test_loss_curr, d_pr = sess.run([G_sample, MSE_train_loss, MSE_test_loss, D_prob],
                        feed_dict = {X: testX, M: testM,
                                     Z: testM * testX + (1-testM) * sample_Z(testM.shape[0], Dim), H: testM, NZ:(testX > 0)*1})
	
        imp = transformBack(preds, testX, testM)
	print(imp.shape, testM.shape, testX.shape, data_full.shape)
	#print([nonZeroMean(imp[:,i]) for i in range(7150, 7155)])
	mse = ((data_full-imp)**2).mean(axis=None)
	mse_preds = ((data_full-np.transpose(preds))**2).mean(axis=None)
	mse_masked = (((data_full-imp)*(1-np.transpose(testM)))**2).mean(axis=None)
	print("threshold:", cutoff_i, cutoffs[cutoff_i], maxIters[cutoff_i], sparsity[cutoff_i])
	print("MSE=", mse)
	nz = (1-testM).sum(axis=1)
	print("AvgNImputsPerCell=", np.min(nz), np.max(nz), np.median(nz))
	print("NZeros=", np.count_nonzero(imp)/(imp.shape[0]*imp.shape[1]*1.0), np.count_nonzero(testX)/(imp.shape[0]*imp.shape[1]*1.0))
	print("NZeros=", np.count_nonzero(1-testM)/(imp.shape[0]*imp.shape[1]*1.0), np.count_nonzero(imp*(1-np.transpose(testM)))*1.0/np.count_nonzero(1-testM))

	percent_nonzero = np.count_nonzero(imp*(1-np.transpose(testM)))*1.0/np.count_nonzero(1-testM)

        #np.savetxt('imputation_gain_data/'+'/gain_preds_pr.txt', d_pr, delimiter=',', fmt="%f")
	
        lib.plot.plot('imputation_gain_data/'+'/loss', "MSE_train", MSE_train_loss_curr)
        lib.plot.plot('imputation_gain_data/'+'/loss', "MSE_test", MSE_test_loss_curr)
	lib.plot.plot('imputation_gain_data/'+'/MSE', "MSE", mse)
	lib.plot.plot('imputation_gain_data/'+'/MSE', "MSE preds", mse_preds)
	lib.plot.plot('imputation_gain_data/'+'/MSE', "MSE imp only", mse_masked)
	lib.plot.plot('imputation_gain_data/'+'/NZeros', "NZeros_imp", np.count_nonzero(imp)/(imp.shape[0]*imp.shape[1]*1.0))
	lib.plot.plot('imputation_gain_data/'+'/NZeros', "NZeros_masked_imp", np.count_nonzero(imp*(1-np.transpose(testM)))*1.0/np.count_nonzero(1-testM))
        lib.plot.flush()

    if it % 5000 == 0 and it > 1:
        #imp_final = transformBackAll(preds, testX, testM, filteredGenes, filteredCells, maxPerCell)
        np.savetxt('imputation_gain_data/'+'/gain_'+out_suffix+'_transformed.csv', imp, delimiter=',', fmt="%f")
        #np.savetxt('imputation_gain_data/'+'/gain_'+out_suffix+"_"+str(it)+'.csv', imp_final, delimiter=',', fmt="%f")



    if percent_nonzero > 0.95 or it > maxIters[cutoff_i]:
	np.savetxt('imputation_gain_data/'+'/gain_'+out_suffix+'_'+str(it)+'_'+str(cutoffs[cutoff_i])+'_transformed.csv', imp, delimiter=',', fmt="%f")
	cutoff_i += 1
        mask, geneAvgs = getMask(data_f, cutoffs[cutoff_i], sparsity[cutoff_i])
        trainM = mask
        testM = trainM[:,:]
        trainX = np.transpose(imp)
        testX = trainX[:,:]
	percent_nonzero = 0
	print("\n=========================\nNew Cutoff : ", cutoffs[cutoff_i])
	print("NZeros=", np.count_nonzero(1-testM)/(imp.shape[0]*imp.shape[1]*1.0), np.count_nonzero(imp*(1-np.transpose(testM)))*1.0/np.count_nonzero(1-testM))

        lib.plot.flush()

    if it == 100000:#1000:
        mask, geneAvgs = getMask(data_f, 1)#0.4)
        trainM = mask
        testM = trainM[:,:]
        trainX = np.transpose(imp)
        testX = trainX[:,:]
    if it == 250000:#25000:
        mask, geneAvgs = getMask(data_f, 0.75)#0.25)
        trainM = mask
        testM = trainM[:,:]
        trainX = np.transpose(imp)
        testX = trainX[:,:]
    if it == 500000:#50000:
        mask, geneAvgs = getMask(data_f, 0.5)#0.2)
        trainM = mask
        testM = trainM[:,:]
        trainX = np.transpose(imp)
        testX = trainX[:,:]


    
    lib.plot.tick()
    #if it % 20000 == 1:
    #    saver.save(sess, "imputation_gain_data/gain_model", global_step=it)

    
