""" Simple implementation of Viterbi algorithm for protein secondary structure
prediction. BE 562 course. """

states = ('Helix','Sheet','Turn ','Other')
start_probability = {'Helix': 0.0, 'Sheet': 0.0,'Turn ':0.0,'Other':1.0}
 
transition_probability = {
 'Helix' : {'Helix': 0.7, 'Sheet': 0.2,'Turn ':0.1,'Other':0.0},
 'Sheet' : {'Helix': 0.3, 'Sheet': 0.6,'Turn ':0.0,'Other':0.1},
 'Turn ' : {'Helix': 0.1, 'Sheet': 0.1,'Turn ':0.1,'Other':0.7},
 'Other' : {'Helix': 0.2, 'Sheet': 0.6,'Turn ':0.0,'Other':0.2}
 }
 
emission_probability = {
 'Helix' : {'M': 0.5, 'L': 0.25,'N':0.05,'E':0.15,'A':0.02,'G':0.03},
 'Sheet' : {'M': 0.05, 'L': 0.10,'N':0.35,'E':0.30,'A':0.15,'G':0.05},
 'Turn ' :    {'M': 0.0, 'L': 0.2,'N':0.1,'E':0.05,'A':0.2,'G':0.45},
 'Other' :   {'M': 0.05, 'L': 0.1,'N':0.1,'E':0.25,'A':0.10,'G':0.35}
 }

def print_dptable(V):
  s = "    " + " ".join(("%7d" % i) for i in range(1,len(V)+1)) + "\n"
  for y in V[0]:
      s += "%.5s: " % y
      s += " ".join("%.7s" % ("%f" % v[y]) for v in V)
      s += "\n"
  print(s)

def viterbi(obs, states, start_p, trans_p, emit_p):
  V = [{}]
  path = {}

  # Initialize V and path with start probabilities
  for state in states:
    V[0][state] = start_p[state] * emit_p[state][obs[0]]
    path[state] = [state]

  # Run Viterbi for t > 0
  for t in range(1, len(obs)):
    V.append({})
    newpath = {}

    for state in states:
        #calculate and record maximum probability and corresponding state
        (maxprob, maxstate) = max((V[t-1][tempstate] * trans_p[tempstate][state] * emit_p[state][obs[t]], tempstate) for tempstate in states)
        V[t][state] = maxprob
        newpath[state] = path[maxstate] + [state]

    path = newpath # Reassign current path

  print_dptable(V)

  #return path with highest probability
  (prob, state) = max((V[t][state], state) for state in states)
  return ('Probability:',prob,'States:',path[state])

def forward(obs, states, start_p, trans_p, emit_p):
  f = [{}]
  P = 0

  #Initialization
  for state in states:
    f[0][state] = start_p[state] * emit_p[state][obs[0]]

  #Iteration
  for t in range(1, len(obs)):
    f.append({})

    for state in states:
      probs = []
      for tempstate in states:
        probs.append(f[t-1][tempstate] * trans_p[tempstate][state] * emit_p[state][obs[t]])
      f[t][state] = sum(probs)
    
    #Termination
    for state in states:
      P = P+f[t][state]

  return ('Total P'+str(obs)+' = '+str(P))

observations = ('M', 'L', 'A','E')
print viterbi(observations,states,start_probability,transition_probability,emission_probability)
print forward(observations,states,start_probability,transition_probability,emission_probability)