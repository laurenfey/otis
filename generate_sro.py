import numpy as np
import numpy.random as rand
import sys
import time

print("-----------------------------------------")
print("          ____  _______________")
print("         / __ \\/_  __/  _/ ___/")
print("        / / / / / /  / / \\__ \\")
print("       / /_/ / / / _/ / ___/ /")
print("       \\____/ /_/ /___//____/ ")
print("    Order Through Informed Swapping")

print("                            __")
print("     ,                    ,\" e`--o")
print("    ((                   (  | __,\'")
print("     \\\\~----------------' \_;\/")
print("     (                      /")
print("     /) ._______________.  )")
print("    (( (               (( (")
print("     ``-\'               ``-\')")
print("-----------------------------------------")

#read in input file
filename=sys.argv[1]
try:
    f = open(filename,'r')
except IOError:
    print('Cannot open input file %s',filename)
print('Reading input file %s'%filename)
lattice_type = f.readline().split()[0]
if lattice_type != 'BCC' and lattice_type != 'FCC':
    raise ValueError('lattice type must be BCC or FCC')
else:
    print('Using %s primitive lattice'%lattice_type)
l = f.readline().split()
if len(l) < 3:
    raise ValueError('Must provide three values corresponding to N1 N2 and N3')
N1,N2,N3 = [int(l[i]) for i in range(3)]
if min([N1,N2,N3]) <= 0:
    raise ValueError('N1 N2 N3 must be greater than zero')
else:
    print('N1=%d N2=%d N3=%d\n%d total atoms'%(N1,N2,N3,N1*N2*N3))
lattice_parameter = float(f.readline().split()[0])
components = int(f.readline().split()[0])
if components <= 0:
    raise ValueError('Number of components must be greater than zero')
labels = f.readline().split()
if len(l) < components:
    raise ValueError('Number of element names does not match number of components')
labels = labels[:components]
l = f.readline().split()
if len(l) < components:
    raise ValueError('Concentration values do not match number of components')
conc = np.array(l[:components],dtype=float)
if abs(np.sum(conc) - 1) > 0.1:
    raise ValueError('Sum of concentrations must be 1')
conc /= np.sum(conc)
print('\nConcentrations:')
for i in range(components):
    print('%s: %0.5f '%(labels[i],conc[i]))

prob = np.zeros((components,components))
input_type = f.readline().split()[0]
if input_type == 'wc':
    for i in range(components):
        l = f.readline().split()
        if len(l) < components:
            raise ValueError('Too few values on Warren-Cowley input')
        for j in range(components):
            prob[i,j] = float(l[j])*(int(i==j) - conc[j])+conc[j]
elif input_type == 'prob':
    for i in range(components):
        l = f.readline().split()
        if len(l) < components:
            raise ValueError('Too few values on probability input')
        for j in range(components):
            prob[i,j] = float(l[j])
else:
    raise ValueError('Input type line must be wc or prob')
print('\nGoal Probability Matrix:')
print('\t'+'\t'.join(labels))
for i in range(components):
    print(labels[i] + '\t' + '\t'.join('%0.4f'%prob[i,j] for j in range(components)))
output_file = f.readline().split()[0]
f.close()

#create lattice and neighbor lists
coords = np.zeros((N1*N2*N3,3))
pos = np.zeros((N1*N2*N3,3))
if lattice_type == 'BCC':
    coord_num = 8
    a1 = np.array([1,1,-1])/np.sqrt(3)
    a2 = np.array([-1,1,1])/np.sqrt(3)
    a3 = np.array([1,-1,1])/np.sqrt(3)
else:
    coord_num = 12
    a1 = np.array([1,1,0])/np.sqrt(2)
    a2 = np.array([0,1,1])/np.sqrt(2)
    a3 = np.array([1,0,1])/np.sqrt(2)
neighbors = np.zeros((N1*N2*N3,coord_num),dtype=int)
index = np.zeros((N1*N2*N3),dtype=int)
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i
            coords[idx] = i,j,k
            pos[idx] = i*a1 +  j*a2 + k*a3
            neighbors[idx,0] = N1*N2*k + N1*j + (i+1)%N1
            neighbors[idx,1] = N1*N2*k + N1*j + (i-1)%N1
            neighbors[idx,2] = N1*N2*k + N1*((j+1)%N2) + i
            neighbors[idx,3] = N1*N2*k + N1*((j-1)%N2) + i
            neighbors[idx,4] = N1*N2*((k+1)%N3) + N1*j + i
            neighbors[idx,5] = N1*N2*((k-1)%N3) + N1*j + i
            if lattice_type == 'BCC':
                neighbors[idx,6] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neighbors[idx,7] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + (i-1)%N1
            else:
                neighbors[idx,6] = N1*N2*((k+1)%N3) + N1*j + (i-1)%N1
                neighbors[idx,7] = N1*N2*((k-1)%N3) + N1*j + (i+1)%N1
                neighbors[idx,8] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + i
                neighbors[idx,9] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + i
                neighbors[idx,10] = N1*N2*k + N1*((j-1)%N2) + (i+1)%N1
                neighbors[idx,11] = N1*N2*k + N1*((j+1)%N2) + (i-1)%N1
print('\nlattice and neighbor lists built')

#assign initial element types
neigh_types = np.zeros((len(coords),components),dtype=int)
lattice = np.zeros(len(coords),dtype=int)
for i in range(components):
    a = int(np.sum(conc[:i])*N1*N2*N3)
    b = int(np.sum(conc[:i+1])*N1*N2*N3)
    lattice[a:b] = i
rand.shuffle(lattice)
print('assigned initial atom types')

#set goal and current bond numbers
goal_bonds = np.zeros((components,components))
curr_bonds = np.zeros((components,components))
for i in range(components):
    for j in range(components):
        goal_bonds[i,j] = N1*N2*N3*coord_num*prob[i,j]*conc[i]
print("goal bonds set")

for i in range(len(coords)):
    for j in neighbors[i]:
        curr_bonds[lattice[i],lattice[j]] += 1
        neigh_types[i,lattice[j]] += 1

delta_swap = np.zeros((len(coords),components,components,components))
for i in range(components):
    for j in range(components):
        delta_swap[lattice==i,j,i,:] -= neigh_types[lattice==i]
        delta_swap[lattice==i,j,:,i] -= neigh_types[lattice==i]
        delta_swap[lattice==i,j,j,:] += neigh_types[lattice==i]
        delta_swap[lattice==i,j,:,j] += neigh_types[lattice==i]
print('calculated current bond numbers')

pairs = [(i,j) for i in range(0,components) for j in range(0,components) if i != j]
rand.shuffle(pairs)
swaps = 0
rejected = 0
tol = 0.0001*N1*N2*N3*coord_num/components
current = [[] for i in range(len(pairs))]
init_time = time.time()
print('bond number tolerance = %f'%tol)

print('Swapping atoms')
max_diff = 0
for p in pairs:
    diff = abs(goal_bonds[p] - curr_bonds[p])
    if diff > max_diff:
        max_diff = diff
init_diff = max_diff
progress = 0
while max_diff > tol:
    for ip,p in enumerate(pairs):
        first = np.argwhere((np.sum(np.abs(curr_bonds+delta_swap[:,p[1]]-goal_bonds),axis=(1,2)) < np.sum(np.abs(curr_bonds-goal_bonds))) & (lattice==p[0]))
        if len(first) == 0:
            idx_0 = rand.randint(0,len(coords))
            while lattice[idx_0] != p[0]:
                idx_0 = rand.randint(0,len(coords))
        else:
            idx_0 = first[rand.randint(0,len(first))][0]

        second = np.argwhere((np.sum(np.abs(curr_bonds+delta_swap[:,p[0]]+delta_swap[idx_0,p[1]]-goal_bonds),axis=(1,2)) < np.sum(np.abs(curr_bonds-goal_bonds))) & (lattice==p[1]))
        if len(second) == 0:
            idx_1 = rand.randint(0,len(coords))
            while lattice[idx_1] != p[1]:
                idx_1 = rand.randint(0,len(coords))
        else:
            idx_1 = second[rand.randint(0,len(second))][0]

        if np.sum(abs(curr_bonds+delta_swap[idx_0,p[1]]+delta_swap[idx_1,p[0]]-goal_bonds)) <= np.sum(abs(curr_bonds-goal_bonds)):
            swaps += 1
            neighbors_0 = neighbors[idx_0]
            neighbors_1 = neighbors[idx_1]

            for n in neighbors_0:
                neigh_types[n,p[0]] -= 1
                neigh_types[n,p[1]] += 1
                curr_bonds[p[0],lattice[n]] -= 1
                curr_bonds[p[1],lattice[n]] += 1
                curr_bonds[lattice[n],p[0]] -= 1
                curr_bonds[lattice[n],p[1]] += 1
                for j in range(components):
                    delta_swap[n,j,lattice[n],p[1]] -= 1
                    delta_swap[n,j,p[1],lattice[n]] -= 1
                    delta_swap[n,j,lattice[n],p[0]] += 1
                    delta_swap[n,j,p[0],lattice[n]] += 1
                    delta_swap[n,j,j,p[1]] += 1
                    delta_swap[n,j,p[1],j] += 1
                    delta_swap[n,j,j,p[0]] -= 1
                    delta_swap[n,j,p[0],j] -= 1
            for n in neighbors_1:
                neigh_types[n,p[1]] -= 1
                neigh_types[n,p[0]] += 1
                curr_bonds[p[1],lattice[n]] -= 1
                curr_bonds[p[0],lattice[n]] += 1
                curr_bonds[lattice[n],p[1]] -= 1
                curr_bonds[lattice[n],p[0]] += 1
                for j in range(components):
                    delta_swap[n,j,lattice[n],p[1]] += 1
                    delta_swap[n,j,p[1],lattice[n]] += 1
                    delta_swap[n,j,lattice[n],p[0]] -= 1
                    delta_swap[n,j,p[0],lattice[n]] -= 1
                    delta_swap[n,j,j,p[1]] -= 1
                    delta_swap[n,j,p[1],j] -= 1
                    delta_swap[n,j,j,p[0]] += 1
                    delta_swap[n,j,p[0],j] += 1

            if idx_0 in neighbors_1:
                curr_bonds[p] += 2
                curr_bonds[p[1],p[0]] +=2
                curr_bonds[p[0],p[0]] -= 2
                curr_bonds[p[1],p[1]] -= 2

            lattice[idx_0] = p[1]
            lattice[idx_1] = p[0]

            for j in range(components):
                delta_swap[idx_0,j] = np.zeros((components,components))
                delta_swap[idx_0,j,lattice[idx_0],:] -= neigh_types[idx_0]
                delta_swap[idx_0,j,:,lattice[idx_0]] -= neigh_types[idx_0]
                delta_swap[idx_0,j,j,:] += neigh_types[idx_0]
                delta_swap[idx_0,j,:,j] += neigh_types[idx_0]

                delta_swap[idx_1,j] = np.zeros((components,components))
                delta_swap[idx_1,j,lattice[idx_1],:] -= neigh_types[idx_1]
                delta_swap[idx_1,j,:,lattice[idx_1]] -= neigh_types[idx_1]
                delta_swap[idx_1,j,j,:] += neigh_types[idx_1]
                delta_swap[idx_1,j,:,j] += neigh_types[idx_1]

        else:
            rejected +=1

        max_diff = 0
        for pair in pairs:
            diff = abs(goal_bonds[pair[0],pair[1]] - curr_bonds[pair[0],pair[1]])
            if diff > max_diff:
                max_diff = diff

        if swaps % 100 == 0:
            progress = max((init_diff-max_diff)/init_diff,progress)
            sys.stdout.write('\r['+'#'*int(progress*50)+' '*(49-int(progress*50))+'] %d'%(progress*100) +'% complete')
            sys.stdout.flush()


print('\nTime to generate SRO lattice (seconds): %f' % (time.time()-init_time))
print('%d Total Swaps. %0.2f' %(swaps,swaps/(swaps+rejected)*100) + '% Acceptance Rate')

actual_bonds = np.zeros((components,components))
for i in range(len(coords)):
    for j in neighbors[i]:
        actual_bonds[lattice[i],lattice[j]] += 1
for i in range(components):
    actual_bonds[i] /= conc[i]*N1*N2*N3*coord_num

print('\nFinal Probability Matrix:')
print('\t'+'\t'.join(labels))
for i in range(components):
    print(labels[i] + '\t' + '\t'.join('%0.4f'%actual_bonds[i,j] for j in range(components)))

f2 = open(output_file,'w')
f2.write("%d\n"%(len(coords)))
f2.write('Lattice=\"%0.8f %0.8f %0.8f '%tuple(a1*N1*lattice_parameter))
f2.write('%0.8f %0.8f %0.8f '%tuple(a2*N2*lattice_parameter))
f2.write('%0.8f %0.8f %0.8f\" pbc="T T T"\n'%tuple(a3*N3*lattice_parameter))
for i in range(len(coords)):
    f2.write("%s %f %f %f\n"%(labels[lattice[i]],pos[i,0]*lattice_parameter,pos[i,1]*lattice_parameter,pos[i,2]*lattice_parameter))
f2.close()

print("\nLattice written to %s"%output_file)
