import sys
import pandas
if len(sys.argv) != 2:
	print("Usage:", sys.argv[0] + "  input.csv")
else:
	data = pandas.read_csv(sys.argv[1])
	data = data.sort_values('Score', ascending=False)
	data.to_csv(sys.argv[1])
