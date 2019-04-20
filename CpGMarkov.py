#Alan Pacheco
#Generate random sequences from M1, M2 and calculate probabilities of coming from M1, M2.

import random

states = ('A','C','G','T')
start_probability = {'A': 0.25, 'C': 0.25,'G':0.25,'T':0.25}
 
chain_length = 6

M1_transition_probability = {
 'A' : {'A': 0.180, 'C': 0.274,'G':0.426,'T':0.120},
 'C' : {'A': 0.171, 'C': 0.367,'G':0.274,'T':0.188},
 'G' : {'A': 0.161, 'C': 0.339,'G':0.375,'T':0.125},
 'T' : {'A': 0.079, 'C': 0.355,'G':0.384,'T':0.182}
 }

M2_transition_probability = {
 'A' : {'A': 0.300, 'C': 0.205,'G':0.285,'T':0.210},
 'C' : {'A': 0.322, 'C': 0.298,'G':0.078,'T':0.302},
 'G' : {'A': 0.248, 'C': 0.246,'G':0.298,'T':0.208},
 'T' : {'A': 0.177, 'C': 0.239,'G':0.292,'T':0.292}
 }

emission_probability = {
 'A' : {'A': 1, 'C': 0,'G':0,'T':0},
 'C' : {'A': 0, 'C': 1,'G':0,'T':0},
 'G' : {'A': 0, 'C': 0,'G':1,'T':0},
 'T' : {'A': 0, 'C': 0,'G':0,'T':1}
 }

#use the weights in a probability matrix to select a nucleotide
def SelectNuc(d):
    r = random.random()
    s = 0.0
    for k, w in d.iteritems():
        s += w
        if r < s: #if the random number is less than the next cumulative probability
        	return k

#Make the Markov chain
def makechain(states,starts,p_trans):
	chain = []
	for i in range(chain_length):
		if i == 0:
			chain.append(SelectNuc(starts))
		else:
			chain.append(SelectNuc(p_trans[chain[i-1]]))
	return chain

#Calculate the probability of the chain coming from one dist or another
def calc_prob(chain, trans_p, em_p):
	prob = 1
	for i in range(len(chain)-1):
		letter = chain[i]
		nextletter = chain[i+1]

		prob *= trans_p[letter][nextletter] * em_p[letter][letter]
	return prob

M1_or_m2 = []
trials = 10000
for i in range(trials):
	chain = makechain(states,start_probability,M1_transition_probability)
	M1_prob = calc_prob(chain,M1_transition_probability,emission_probability)
	M2_prob = calc_prob(chain,M2_transition_probability,emission_probability)
	if M1_prob < M2_prob:
		M1_or_m2.append(1)
	else:
		M1_or_m2.append(0)

print 'Part a: M2 probability is higher than M1 probability',sum(M1_or_m2),'out of',trials,'times.'

M1_or_m2 = []
trials = 10000
for i in range(trials):
	chain = makechain(states,start_probability,M2_transition_probability)
	M1_prob = calc_prob(chain,M1_transition_probability,emission_probability)
	M2_prob = calc_prob(chain,M2_transition_probability,emission_probability)
	if M1_prob > M2_prob:
		M1_or_m2.append(1)
	else:
		M1_or_m2.append(0)

print 'Part b: M1 probability is higher than M2 probability',sum(M1_or_m2),'out of',trials,'times.'