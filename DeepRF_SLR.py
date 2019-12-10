#
# main script to perform deep reinforcement learning designed spin-echo pulse design (DeepRF_SLR)
# implemented by Dongmyung Shin, Seoul National University
# shinsae11@gmail.com, http://list.snu.ac.kr
#
# install all pre-requirements before running this script
# follow this URL for info.: http://github.com/SNU-LIST/DeepRF_SLR
#
# usage: e.g., type this command in Ubuntu terminal:
# python3 DeepRF_SLR.py 3 6 6 --gpu 0 --lr 1e-4 --nn 8 --nl 256 --gamma 1.0 --iter 100000
#

#%% import libraries
import os
import matlab.engine
import numpy as np
import tensorflow as tf
from scipy.io import savemat
import argparse
import time

#%% arguments
parser = argparse.ArgumentParser()
parser.add_argument('NB', type = int, help = 'Number of bands')
parser.add_argument('TBW', type = int, help = 'Time-bandiwdth product')
parser.add_argument('BS', type = int, help = 'N times of slice thickenss (e.g. 6->500%)')
parser.add_argument('--gpu', type = int, default = 0, help = 'activiated GPU number when having multiple GPUs (default=0)')
parser.add_argument('--lr', type = float, default = 1e-4, help = 'learning rate (default=1e-4)')
parser.add_argument('--nn', type = int, default = 256, help = 'number of node for each hidden layer (default=256)')
parser.add_argument('--nl', type = int, default = 8, help = 'number of layers except input layer and including output layer (default=8) (>2)')
parser.add_argument('--gamma', type = float, default = 1.0, help = 'discount rate (1.0) (<=1.0)')
parser.add_argument('--iter', type = int, default = 100000, help = 'number of flipping until termination (default=100000)')
args = parser.parse_args()

#%% pulse design parameters
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
n = 512.0;                          # number of time points
nb = float(args.NB);                # number of bands
tb = float(args.TBW);               # Time-bandwidth product
bandsep = float(args.BS)*tb;        # Slice gap is BS x slice thickness
d1e = 0.01;                         # combined Mxy ripple, passband
d2e = 0.01;                         # combined Mxy ripple, stopband

#%% initial minimum phase SLR pulse design
eng = matlab.engine.start_matlab()
bp,d1,tbrf,r,N = eng.min_phase_design(n,tb,d1e,d2e,bandsep,nb,nargout=5)
state_size,idxPass,r = eng.num_root(bp,tbrf,r,N,nargout=3)
state_size = int(state_size);
action_size = state_size
print("number of passband roots: {}".format(state_size))

#%% learning hyperparameters
learning_rate = float(args.lr)              # learning rate
num_node = int(args.nn)                     # number of node in eacy layer
num_layer = int(args.nl)                    # number of layers
gamma = float(args.gamma)                   # discount factor
max_steps = action_size                     # number of steps in one episode
max_iter = int(args.iter)                   # termination condition

#%% function definitions
def discount_and_nomarlize_rewards(episode_rewards):
    discounted_episode_rewards = np.zeros_like(episode_rewards)
    cumulative = 0.0
    for i in reversed(range(len(episode_rewards))):
        cumulative = cumulative * gamma + episode_rewards[i]
        discounted_episode_rewards[i] = cumulative
    return discounted_episode_rewards

def fc_layer_with_leaky(input, channels_in, channels_out, name = "fc"):
    w = tf.Variable(tf.truncated_normal([channels_in, channels_out], stddev = \
                                        tf.sqrt( tf.divide(2.0, channels_in) ))) # He initializer
    b = tf.Variable(tf.constant(0.0, shape = [channels_out]))
    act = tf.nn.leaky_relu(tf.matmul(input, w) + b)
    return act
    
def fc_layer_without_act(input, channels_in, channels_out, name = "fc"):
    w = tf.Variable(tf.truncated_normal([channels_in, channels_out], stddev = \
                                        tf.sqrt( tf.divide(2.0, channels_in) ))) # He initializer
    b = tf.Variable(tf.constant(0.0, shape = [channels_out]))
    act = tf.matmul(input, w) + b
    return act    

#%% deep neural network definition
input_ = tf.placeholder(tf.float32,[None,int(n)])
actions = tf.placeholder(tf.float32,[None,action_size])
advantages_ = tf.placeholder(tf.float32,[None,])
fc = fc_layer_with_leaky(input_, int(n), num_node)
for layer in range(num_layer - 2):
    fc = fc_layer_with_leaky(fc, num_node, num_node)
fc = fc_layer_without_act(fc, num_node, action_size)
action_distribution = tf.nn.softmax(fc) + 1e-8
log_prob = tf.reduce_sum(tf.multiply(actions,tf.log(action_distribution)), [1])
loss = -tf.reduce_mean(tf.multiply(log_prob,advantages_))
train_opt = tf.train.AdamOptimizer(learning_rate).minimize(loss)

#%% run policy gradient method
reward_list = []
loss_list = []
val_loss_list = []
rf_list = []
amp_list = []
iter_list = []
time_list = []
pattern_list = []

# track the best |RF| for all iterations
bestRF_global = 1
iter = 0

# arrays for saving the results
episode_states,episode_actions,episode_action_,episode_rewards = [],[],[],[]
episode_states_rf,episode_pattern,episode_val = [],[],[]

# non-flipped RF design
initial_rf_mag,initial_rf_angle,_ = eng.update_pulse(r,matlab.double(np.zeros((state_size,1)).tolist()),idxPass,N,d1,nargout=3)
initial_rf = np.array(initial_rf_mag)*np.exp(1j*np.array(initial_rf_angle))
initial_rf_concat = np.real(initial_rf)
initial_rf_max = np.max(np.abs(initial_rf))
rf_list.append(initial_rf_concat)
amp_list.append(initial_rf_max)
iter_list.append(0)

# start clock
start_time = time.time()

with tf.Session() as sess:
            
    sess.run(tf.global_variables_initializer())
    
    for episode in range(10000000000):
        
        # Launch the game
        state = np.zeros((1,state_size))
        state_rf = np.copy(initial_rf_concat)
        
        # track the best |RF| in one episode for policy gradient update
        bestRF_in_episode = 1
        bestState_in_episode = np.copy(state)
        bestStateRF_in_episode = np.copy(state_rf)
        
        for iteration in range(max_steps + 1):

            # start of greedy tree search
            if iteration == max_steps:
                      
                pattern = np.copy(bestState_in_episode)
                pattern_rf = np.copy(bestStateRF_in_episode)
                                                
                prevRF = bestRF_in_episode
                print("start |RF|: ",prevRF)
                
                while 1:
                                            
                    max_root = 0
                    max_r = 0
                    
                    # span all actions
                    for a_num in range(state_size):
                        temp_pattern = np.copy(pattern)
                        temp_pattern[0,a_num] = abs(temp_pattern[0,a_num]-1)
                        rfmag,rfangle,rfroots = eng.update_pulse(r,matlab.double(temp_pattern.tolist()),idxPass,N,d1,nargout=3)
                        state_rf_temp = np.array(rfmag)*np.exp(1j*np.array(rfangle))
                        temp_pattern_rf = np.concatenate((np.real(state_rf_temp),np.imag(state_rf_temp)),axis=1)
                        currentRF = np.max(rfmag)
                        reward = prevRF - currentRF
                        
                        iter += 1
                        
                        if reward > max_r:
                            max_r = reward
                            max_root = a_num                                                                           
                            pattern_ = np.copy(temp_pattern)
                            pattern_rf_ = np.copy(temp_pattern_rf)
                            prevRF_ = np.max(np.abs(state_rf_temp))
                            
                        if bestRF_global > currentRF:
                            bestRF_global = currentRF
                            amp_list.append(bestRF_global)
                            iter_list.append(iter)
                            rf_list.append(state_rf_temp)
                            hold_time = time.time()
                            time_list.append(hold_time - start_time)
                            pattern_list.append(temp_pattern)
                        
                    if max_r == 0:
                        break
                    else:
                        pattern = np.copy(pattern_)
                        pattern_rf = np.copy(pattern_rf_)
                        prevRF = prevRF_
                        print("updated |RF|: ",prevRF_)
            
            else:                
            
                action_prob_dist = sess.run(action_distribution,feed_dict={input_:state_rf})
                
                # enforce flipped action do not selected
                action_prob_dist[state == 1] = 0.0
                action_prob_sum_to_one = action_prob_dist / action_prob_dist.sum(axis=1)
                
                action = np.random.choice(range(action_prob_sum_to_one.shape[1]),p=action_prob_sum_to_one.ravel())
                iter += 1
                
                episode_states_rf.append(state_rf[0,:])
                
                # State transition
                state[0,action] = abs(state[0,action]-1)
                rfmag,rfangle,_ = eng.update_pulse(r,matlab.double(state.tolist()),idxPass,N,d1,nargout=3)
                state_rf_temp = np.array(rfmag)*np.exp(1j*np.array(rfangle))
                state_rf = np.real(state_rf_temp)
                maxRF = np.max(np.abs(state_rf_temp))                
                
                # if new minimum peak is found
                if bestRF_global > maxRF:
                    bestRF_global = maxRF
                    amp_list.append(bestRF_global)
                    rf_list.append(state_rf_temp)
                    iter_list.append(iter)                    
                    hold_time = time.time()
                    time_list.append(hold_time - start_time)
                    pattern_list.append(state)
                    
                if bestRF_in_episode > maxRF:
                    bestRF_in_episode = maxRF
                    bestState_in_episode = np.copy(state)
                    bestStateRF_in_episode = np.copy(state_rf)                        
                    
                action_ = np.zeros(action_size)
                action_[action] = 1
                
                episode_actions.append(action_)
                episode_action_.append(action)
                
                if iteration == max_steps - 1:
                    reward = 1/bestRF_in_episode - 1/initial_rf_max
                    episode_rewards.append(reward)                        
                else:
                    episode_rewards.append(0)                    
                    
        # track the reward
        reward_list.append(np.sum(episode_rewards)) # equal to best result in one episode
        
        # Calculate discounted reward
        discounted_episode_rewards = discount_and_nomarlize_rewards(episode_rewards) # no discount            
        advantages = discounted_episode_rewards

        # train network
        log_prob_temp,loss_temp,_ = sess.run([log_prob,loss,train_opt],feed_dict={
                           input_:np.vstack(np.array(episode_states_rf)),
                           actions:np.vstack(np.array(episode_actions)),
                           advantages_:advantages})                
        loss_list.append(loss_temp)
        
        print("Episode: ",episode+1)
        print("Iterations: ",iter+1)
        print("Best max |RF|: ",bestRF_global)
        print("Reward Sum: ",np.sum(episode_rewards))
        print("Network loss: ",loss_temp)
        print("===============================================================")
        
        if iter > max_iter:  
            hold_time = time.time()
            time_list.append(hold_time - start_time)
            savemat('DeepRF_SLR_refo_design.mat',
                                                dict([('reward_list',reward_list),
                                                      ('loss_list',loss_list),
                                                      ('rf_list',rf_list),
                                                      ('amp_list',amp_list),
                                                      ('iter_list',iter_list),
                                                      ('tbrf',tbrf),
                                                      ('state_size',state_size),
                                                      ('action_size',action_size),
                                                      ('time_list',time_list),
                                                      ('pattern_list',pattern_list),
                                                      ('args',args)]))
            break
        
        # Reset the list        
        episode_states,episode_actions,episode_rewards,episode_action_ = [],[],[],[]
        episode_states_rf,episode_pattern,episode_val = [],[],[]
    
