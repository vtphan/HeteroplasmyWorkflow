import csv
import locale

def get_entries(csv_file):
	locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' ) 
	entries = []
	with open(csv_file) as file:
		reader = csv.DictReader(file)
		for row in reader:
			if row['Type'] in ['rRNA', 'tRNA', 'CDS']:
				e = (locale.atoi(row['Minimum']), locale.atoi(row['Maximum']), row['Name'], row['Type'], row['# Intervals'], row['Direction'])
				entries.append(e)
	return sorted(entries)

def search(entries, pos):
	def search_r(L, R, pos):
		if L > R:
			return -1
		m = (L+R)//2
		cur = entries[m]
		if cur[0] <= pos <= cur[1]:
			return m
		if pos < cur[0]:
			return search_r(L, m-1, pos)
		else:
			return search_r(m+1, R, pos)
	return search_r(0, len(entries)-1, pos)

if __name__ == '__main__':
	import sys
	if len(sys.argv) != 2:
		print("Usage:", sys.argv[0] + "  input.csv")
		sys.exit(1)
	E = get_entries(sys.argv[1])
	P = [0,3,50,200,300,1342,1325,9027,9049,154292,155000,155767,160000]
	for p in P:
		pos = search(E,p)
		if pos==-1:
			print(p, "Not found")
		else:
			print(p, E[pos])