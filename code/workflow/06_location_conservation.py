import pandas
import sys
import math

MIN_SCORE = 0.05

#------------------------------------------------------------------------------
def cosine_sim(a,b):
	dot = sum([ a[i]*b[i] for i in range(len(a)) ])
	norm_a = math.sqrt( sum( [a[i]*a[i] for i in range(len(a))] ) )
	norm_b = math.sqrt( sum( [b[i]*b[i] for i in range(len(b))] ) )
	return dot / (norm_a * norm_b)

#------------------------------------------------------------------------------
# https://en.wikipedia.org/wiki/Hellinger_distance
#------------------------------------------------------------------------------
def hellinger_sim(a,b):
	# 1/sqrt(2)
	od2 = 0.7071067811865475
	d = od2 * math.sqrt(sum( (math.sqrt(a[i])-math.sqrt(b[i]))**2 for i in range(len(a))))
	return 1.0 - d

#------------------------------------------------------------------------------
def group_conservation(rows, sim_func):
	if len(rows) <= 1:
		return MIN_SCORE
	sim = 0
	for i in range(len(rows)):
		for j in range(i+1, len(rows)):
			a = (rows[i]['A'], rows[i]['C'], rows[i]['G'], rows[i]['T'], rows[i]['D'])
			b = (rows[j]['A'], rows[j]['C'], rows[j]['G'], rows[j]['T'], rows[j]['D'])
			sim += sim_func(a,b)
	return 2.0 * sim / (len(rows)*(len(rows)-1))

#------------------------------------------------------------------------------
def location_conservation(df, dist_func):
	g = df[['Coordinate','A','C','G','T','D']].groupby('Coordinate')
	score = {}
	for gid in g.groups:
		group = g.get_group(gid)   # group is a df
		rows = [ r[1] for r in group.iterrows() ]
		score[gid] = group_conservation(rows, dist_func)
	return score

#------------------------------------------------------------------------------
if len(sys.argv) != 3:
	print('\n\tUsage: python', sys.argv[0], 'heteroplasmy_file.csv dist')
	print('\tdist is either cosine or hellinger')
	print('\tOutput: conservation scores of heteroplasmic sites\n')
	sys.exit(0)

df = pandas.read_csv(sys.argv[1])
if sys.argv[2] == 'cosine':
	sim_func = cosine_sim
elif sys.argv[2] == 'hellinger':
	sim_func = hellinger_sim
else:
	print('Error: unknown distance function (%s)' % sys.argv[2])
	sys.exit(0)

scores = location_conservation(df, sim_func)
print('Coordinate,Score')
for v in sorted(scores.items()):
	print('%s,%s' % (v[0],v[1]))
