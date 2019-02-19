import csv
import os
import sys

#------------------------------------------------------------------------------
def get_csvfiles(dir):
	return [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') ]

#------------------------------------------------------------------------------
def get_csvfiles_by_nameList(dir, namelist):

	def getFile(f_dir, name):
		for i in f_dir:
			if name in i:
				return i

		return None

	f = []
	f_dir = [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') ]
	with open(namelist,'r') as nl:
		for line in nl:
			if "," in line:
				name, organism = line.strip().split(",")
			else:
				name = line.strip()

			if getFile(f_dir, name):
					f.append(getFile(f_dir, name))

	return f

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

			a, c, g, t, d, i, total = float(row['A']), float(row['C']), float(row['G']), float(row['T']), float(row['D']), float(row['I']), float(row['Total'])
			percentages = sorted([a/total, c/total, g/total, t/total, d/total, i/total])
			highest, second_highest = percentages[5], percentages[4]
			
			if score >= score_threshold:
				high_score.append((pos, score, second_highest, row['GeneProduct'], a/total, c/total, g/total, t/total, d/total, i/total, a, c, g, t, d, i, total))
				# print("%d\t%.4f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f" % (pos, score, a,c,g,t,total,percentages[0],percentages[1],percentages[2],percentages[3]))
			else:
				break
	selected = [ x for x in high_score if x[2] >= percentage_threshold ]
	# ret_val = { s[0] : person_id for s in selected }, { s[0] : s for s in selected }
	return { s[0] : s for s in selected }

#------------------------------------------------------------------------------
def intersect(csvfiles, score_threshold, percentage_threshold, output_file):
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
	scatter_plot([get_individual_id(f) for f in csvfiles], positions, output_file)

#------------------------------------------------------------------------------
def scatter_plot(ids, positions, output_file):
	points = []
	with open(output_file, 'w') as f:
		# print('Coordinate,Sample,Name,GP,A,C,G,T,D,I,CountA,CountC,CountG,CountT,CountD,CountI,total,d')
		f.write('Coordinate,Sample,Name,GP,A,C,G,T,D,I,CountA,CountC,CountG,CountT,CountD,CountI,total,d\n')

		for i,cur_id in enumerate(ids):
			items = []
			for pos,profile in positions.items():
				if cur_id in profile[0]:
					idx = profile[0].index(cur_id)
					item = profile[1][idx]
					# print('%d,%d,%s' % (pos,i+1,cur_id))
					# print(item)
					# print(profile)
					items.append([pos,i+1,cur_id,item[3],item[4],item[5],item[6],item[7],item[8],item[9],item[10],item[11],item[12], item[13], item[14], item[15],item[16],0])

			# Compute the distance to the nearest neighbors
			items.sort()
			if len(items) > 1:
				items[0][-1] = items[1][0] - items[0][0]
				for j in range(1,len(items)-1):
					items[j][-1] = min(items[j][0]-items[j-1][0],  items[j+1][0]-items[j][0])
				items[-1][-1] = items[-1][0] - items[-2][0]

			# for x in items:
			# 	print('%d,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%d,%d' % tuple(x))

			for x in items:
				f.write('%d,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%d,%d\n' % tuple(x))

	print("Finish selecting sites.\n")

def process(params):
	csv_dir = params['csv_dir']
	score_threshold = int(params['score_threshold'])
	percentage_threshold = float(params['percentage_threshold'])
	name_list = params['name_list']
	organellar_type = params['organellar_type']
	result_dir = params['result_dir']

	if name_list:
		files = get_csvfiles_by_nameList(csv_dir, name_list)
	else:
		files = get_csvfiles(csv_dir)

	if not os.path.exists(result_dir):
		os.makedirs(result_dir)

	output_file = os.path.join(result_dir, organellar_type+"_heteroplasmy.csv")

	intersect(files, score_threshold, percentage_threshold, output_file)

	return output_file
	

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# organellar type is either chloroplast or mitochondria
	if len(sys.argv) != 6 and len(sys.argv) != 7:
		print("USAGE: ", sys.argv[0], "  csv_dir score_threshold percentage_threshold organellar_type result_dir")
		print("or")
		print("USAGE: ", sys.argv[0], "  csv_dir score_threshold percentage_threshold name_list.csv organellar_type result_dir")
		sys.exit(0)

	if len(sys.argv) == 6:
		name_list = None
		organellar_type = sys.argv[4]
		result_dir = sys.argv[5]
	else:
		name_list = sys.argv[4]
		organellar_type = sys.argv[5]
		result_dir = sys.argv[6]

	params = {
		'csv_dir': sys.argv[1],
		'score_threshold': int(sys.argv[2]),
		'percentage_threshold': float(sys.argv[3]),
		'name_list': name_list,
		'organellar_type': organellar_type,
		'result_dir': result_dir
	}

	# if len(sys.argv) == 4:
	# 	files = get_csvfiles(sys.argv[1])
	# else:
	# 	files = get_csvfiles_by_nameList(sys.argv[1], sys.argv[4])
	
	# score_threshold = int(sys.argv[2])
	# percentage_threshold = float(sys.argv[3])
	# # filter_csvfile(files[0], 'SRR2147184')
	# intersect(files, score_threshold, percentage_threshold)
	process(params)
