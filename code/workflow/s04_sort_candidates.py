import sys
import pandas


def process(**params):
	input_csv = params['out_csv']
	data = pandas.read_csv(input_csv)
	data = data.sort_values('Score', ascending=False)
	data.to_csv(input_csv)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("Usage:", sys.argv[0] + "  input.csv")
	else:
		data = pandas.read_csv(sys.argv[1])
		data = data.sort_values('Score', ascending=False)
		data.to_csv(sys.argv[1])
