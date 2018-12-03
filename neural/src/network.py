"""
network.py
~~~~~~~~~~

A module to implement the stochastic gradient descent learning
algorithm for a feedforward neural network.  Gradients are calculated
using backpropagation.  Note that I have focused on making the code
simple, easily readable, and easily modifiable.  It is not optimized,
and omits many desirable features.
"""

import pdb
#### Libraries
# Standard library
import random

# Third-party libraries
import numpy as np

class Network(object):

	def __init__(self, sizes):
		"""The list ``sizes`` contains the number of neurons in the
		respective layers of the network.  For example, if the list
		was [2, 3, 1] then it would be a three-layer network, with the
		first layer containing 2 neurons, the second layer 3 neurons,
		and the third layer 1 neuron.  The biases and weights for the
		network are initialized randomly, using a Gaussian
		distribution with mean 0, and variance 1.  Note that the first
		layer is assumed to be an input layer, and by convention we
		won't set any biases for those neurons, since biases are only
		ever used in computing the outputs from later layers."""
		self.num_layers = len(sizes)
		self.sizes = sizes
		self.biases = [np.random.randn(y, 1) for y in sizes[1:]]
		self.weights = [np.random.randn(y, x)
						for x, y in zip(sizes[:-1], sizes[1:])]

	def feedforward(self, a):
		"""Return the output of the network if ``a`` is input."""
		for b, w in zip(self.biases, self.weights):
			a = sigmoid(np.dot(w, a)+b)
		return a

	def SGD(self, training_data, epochs, mini_batch_size, eta,
			test_data=None):
		"""Train the neural network using mini-batch stochastic
		gradient descent.  The ``training_data`` is a list of tuples
		``(x, y)`` representing the training inputs and the desired
		outputs.  The other non-optional parameters are
		self-explanatory.  If ``test_data`` is provided then the
		network will be evaluated against the test data after each
		epoch, and partial progress printed out.  This is useful for
		tracking progress, but slows things down substantially."""
		if test_data: n_test = len(test_data)
		categorized_training_data = {}
		for t in training_data:
			key = int(t[-1])
			if key not in categorized_training_data.keys():
				categorized_training_data[key] = []
			categorized_training_data[key].append(t)
		keys = sorted(categorized_training_data.keys())
		categorized_training_data_lengths = [len(categorized_training_data[k]) for k in keys]
		for j in xrange(epochs):
			mini_batch = [[] for i in keys]
			mini_batch_deltas = [[] for i in keys]
			mini_batch_deltas_s = [0 for i in keys]
			mini_batch_deltas_s2 = [0 for i in keys]
			mini_batch_variances = [0 for i in keys]
			pdb.set_trace()
			for k in keys:
				mini_batch[k].append(random.choice(categorized_training_data[k]))
				mini_batch[k].append(random.choice(categorized_training_data[k]))
				#RAWR
				mini_batch_deltas[k] = [self.backprop(x, y) for x, y, z in mini_batch[k]]
				mini_batch_deltas_s[k] = (sum([a[0] for a,b in mini_batch_deltas[k]]),sum([a[1] for a,b in mini_batch_deltas[k]]))
				mini_batch_deltas_s2[k] = (sum([kk[0]**2 for kk in mini_batch_deltas[k]]),sum([kk[0]**2 for kk in mini_batch_deltas[k]]))
				mini_batch_variances[k] = mini_batch_deltas_s2[k] - mini_batch_deltas_s[k]**2*1.0/2
			
			mini_batch = [random.choice(categorized_training_data[random.randint(0,9)]) for k in range(mini_batch_size)]
			delta_nablas = [self.backprop(x, y) for x, y, z in mini_batch]
			self.update_mini_batch(delta_nablas, eta)
			if test_data:
				print "Epoch {0}: {1} / {2}".format(
					j, self.evaluate(test_data), n_test)
			else:
				print "Epoch {0} complete".format(j)

	def update_mini_batch(self, delta_nablas, eta):
		"""from the delta_nabla_b and delta_nabla_w, Update the network's weights and biases
		 ``eta`` is the learning rate."""
		nabla_b = [np.zeros(b.shape) for b in self.biases]
		nabla_w = [np.zeros(w.shape) for w in self.weights]
		for x, y in delta_nablas:
			nabla_b = [nb+dnb for nb, dnb in zip(nabla_b, x)]
			nabla_w = [nw+dnw for nw, dnw in zip(nabla_w, y)]
		self.weights = [w-(eta/len(delta_nablas))*nw
						for w, nw in zip(self.weights, nabla_w)]
		self.biases = [b-(eta/len(delta_nablas))*nb
					   for b, nb in zip(self.biases, nabla_b)]

	def backprop(self, x, y):
		"""Return a tuple ``(nabla_b, nabla_w)`` representing the
		gradient for the cost function C_x.  ``nabla_b`` and
		``nabla_w`` are layer-by-layer lists of numpy arrays, similar
		to ``self.biases`` and ``self.weights``."""
		nabla_b = [np.zeros(b.shape) for b in self.biases]
		nabla_w = [np.zeros(w.shape) for w in self.weights]
		# feedforward
		activation = x
		activations = [x] # list to store all the activations, layer by layer
		zs = [] # list to store all the z vectors, layer by layer
		for b, w in zip(self.biases, self.weights):
			z = np.dot(w, activation)+b
			zs.append(z)
			activation = sigmoid(z)
			activations.append(activation)
		# backward pass
		delta = self.cost_derivative(activations[-1], y) * \
			sigmoid_prime(zs[-1])
		nabla_b[-1] = delta
		nabla_w[-1] = np.dot(delta, activations[-2].transpose())
		# Note that the variable l in the loop below is used a little
		# differently to the notation in Chapter 2 of the book.  Here,
		# l = 1 means the last layer of neurons, l = 2 is the
		# second-last layer, and so on.  It's a renumbering of the
		# scheme in the book, used here to take advantage of the fact
		# that Python can use negative indices in lists.
		for l in xrange(2, self.num_layers):
			z = zs[-l]
			sp = sigmoid_prime(z)
			delta = np.dot(self.weights[-l+1].transpose(), delta) * sp
			nabla_b[-l] = delta
			nabla_w[-l] = np.dot(delta, activations[-l-1].transpose())
		return (nabla_b, nabla_w)

	def evaluate(self, test_data):
		"""Return the number of test inputs for which the neural
		network outputs the correct result. Note that the neural
		network's output is assumed to be the index of whichever
		neuron in the final layer has the highest activation."""
		test_results = [(np.argmax(self.feedforward(x)), y)
						for (x, y) in test_data]
		return sum(int(x == y) for (x, y) in test_results)

	def cost_derivative(self, output_activations, y):
		"""Return the vector of partial derivatives \partial C_x /
		\partial a for the output activations."""
		return (output_activations-y)

#### Miscellaneous functions
def sigmoid(z):
	"""The sigmoid function."""
	return 1.0/(1.0+np.exp(-z))

def sigmoid_prime(z):
	"""Derivative of the sigmoid function."""
	return sigmoid(z)*(1-sigmoid(z))
