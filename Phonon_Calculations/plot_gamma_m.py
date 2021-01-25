import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

plt.rcParams['axes.labelsize'] = 25
plt.rcParams['axes.titlesize'] = 25
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

natom = 42
nqpt = 11

n_to_plot = 20

fil = 'test_1X1'

f = open(fil,'r')
data = f.readlines()
f.close()

counter = 0
x = [i for i in range(nqpt)]
for j in range(n_to_plot):
    freq = []
    for i in range(nqpt):
        d = data[2+i*(3*natom+1)+j].split()
        freq.append(eval(d[3]))
    plt.plot(x,freq,c='b')
plt.xlabel(' ')
plt.ylabel(r'Frequency (cm$^{-1}$)')

label=(r'$\Gamma $',r'M')
k_nodes = [x[0],x[len(x)-1]]
plt.xticks(k_nodes, label)
#plt.ylim(-2,10)
plt.xlim(0,len(x)-1)
plt.savefig("9d43.png")
#plt.show()
