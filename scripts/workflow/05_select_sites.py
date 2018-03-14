import csv
import os
import sys

#------------------------------------------------------------------------------
def get_csvfiles(dir):
	return [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') ]

#------------------------------------------------------------------------------
def get_individual_id(csvfile):
	filename = csvfile.split('/')[-1]
	prefix = filename.split('.')[0]
	name = prefix.split('_')[0]
	return name

#------------------------------------------------------------------------------
def filter_csvfile(csvfile, score_threshold, percentage_threshold):
	# print("Processing", csvfile)
	high_score = []
	with open(csvfile) as f:
		reader = csv.DictReader(f)
		for row in reader:
			pos, score = int(row['Pos']), float(row['Score'])
			a, c, g, t, total = float(row['A']), float(row['C']), float(row['G']), float(row['T']), float(row['Total'])
			percentages = sorted([a/total, c/total, g/total, t/total])
			highest, second_highest = percentages[3], percentages[2]
			if score >= score_threshold:
				high_score.append((pos, score, second_highest, row['GeneProduct'], a/total, c/total, g/total, t/total, a, c, g, t, total))
				# print("%d\t%.4f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f" % (pos, score, a,c,g,t,total,percentages[0],percentages[1],percentages[2],percentages[3]))
			else:
				break
	selected = [ x for x in high_score if x[2] >= percentage_threshold ]
	# ret_val = { s[0] : person_id for s in selected }, { s[0] : s for s in selected }
	return { s[0] : s for s in selected }

#------------------------------------------------------------------------------
def intersect(csvfiles, score_threshold, percentage_threshold):
	files = {}
	positions = {}
	for f in csvfiles:
		person_id = get_individual_id(f)
		pos = filter_csvfile(f, score_threshold, percentage_threshold)
		for p,profile in pos.items():
			if p not in positions:
				positions[p] = ([], [])
			positions[p][0].append(person_id)
			positions[p][1].append(profile)

	# for k,v in positions.items():
	# 	print(k,v)
	scatter_plot([get_individual_id(f) for f in csvfiles], positions)

#------------------------------------------------------------------------------
def scatter_plot(ids, positions):
	points = []
	print('Coordinate,Sample,Name,GP,A,C,G,T,CountA,CountC,CountG,CountT,total,d')
	for i,cur_id in enumerate(ids):
		items = []
		for pos,profile in positions.items():
			if cur_id in profile[0]:
				idx = profile[0].index(cur_id)
				item = profile[1][idx]
				# print('%d,%d,%s' % (pos,i+1,cur_id))
				# print(profile)
				items.append([pos,i+1,cur_id,item[3],item[4],item[5],item[6],item[7],item[8],item[9],item[10],item[11],item[12],0])

		# Compute the distance to the nearest neighbors
		items.sort()
		if len(items) > 1:
			items[0][-1] = items[1][0] - items[0][0]
			for j in range(1,len(items)-1):
				items[j][-1] = min(items[j][0]-items[j-1][0],  items[j+1][0]-items[j][0])
			items[-1][-1] = items[-1][0] - items[-2][0]

		for x in items:
			print('%d,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%d' % tuple(x))

#------------------------------------------------------------------------------
if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("USAGE: ", sys.argv[0], "  csv_dir score_threshold percentage_threshold")
		sys.exit(0)
	files = get_csvfiles(sys.argv[1])
	# filter_csvfile(files[0], 'SRR2147184')
	score_threshold = float(sys.argv[2])
	percentage_threshold = float(sys.argv[3])
	intersect(files, score_threshold, percentage_threshold)
