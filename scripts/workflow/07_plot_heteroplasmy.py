from bokeh.models import HoverTool, NumeralTickFormatter, FuncTickFormatter, FixedTicker
from bokeh.models import ColumnDataSource, LabelSet, TapTool, Spacer
from bokeh.models import LinearColorMapper, Range1d, Circle
from bokeh.models.widgets import AutocompleteInput, Button, TextInput, Slider
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Blues9, OrRd9, YlOrRd9, Accent8, BuGn9, Set1
from bokeh.layouts import row, column, widgetbox
from bokeh import events

import pandas
import numpy
import sys
import argparse
import csv
import locale
import numpy as np
import math

#------------------------------------------------------------------------------
# Used to set y range of heteroplasmy plot
# Plots (1,0) and (1,1)
#------------------------------------------------------------------------------
VISIBLE_SAMPLE_RANGE = (0,36)
MAX_X = 1

#------------------------------------------------------------------------------
# The entire figure is a 2x3 grid
# DIM[row,column] = (width,height) of plot at location (row,column)
#------------------------------------------------------------------------------
DIM = {
	(0,0) : (1050, 70),
	(0,1) : ( 120, 70),
	(1,0) : (1050,550),
	(1,1) : ( 120,550),
	(2,0) : (1050, 90),
	(2,1) : ( 100, 90),
}

#------------------------------------------------------------------------------
# GENE_INTERVAL[gene_symbol] = [min, max, min_0, max_0, min_1, max_1, ... ]
# Used to label gene products returned by gene zooming search
# Plot (2,0)
#------------------------------------------------------------------------------
GENE_INTERVAL = {}

#------------------------------------------------------------------------------
# HETEROPLASMY_PROBABILITIES[coord] = dict of [(sample_id, probs), .... ]
# Used to plot hbar of prob distributions of heteroplasmy of samples
# Plot (1,1)
#------------------------------------------------------------------------------
HETEROPLASMY_PROBABILITIES = {}

#------------------------------------------------------------------------------
# Command line arguments
#------------------------------------------------------------------------------
ARGS = None

#------------------------------------------------------------------------------
def get_cmd_args():
	parser = argparse.ArgumentParser(description='Create heteroplasmy plot')
	parser.add_argument('genome_name', help='Name of genome (a string)')
	parser.add_argument('genome_annotations', help='Annotations of gene products (csv)')
	parser.add_argument('heteroplasmies', help='Heteroplasmies file (csv)')
	parser.add_argument('conserved_scores', help='conserved scores file (csv)')
	parser.add_argument('output', help='Output file')
	return parser.parse_args()

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------
def main():
	global ARGS
	ARGS = get_cmd_args()

	# Plot heteroplasmic sites
	plasmy_fig, plasmy_source = plot_heteroplasmies()

	annotation_fig = plot_genome_annotations(plasmy_fig)
	label_source, line_source = plot_search_result_annotations(annotation_fig)

	# Plot heteroplasmy probability figure and conservation annotation (in main figure)
	prob_fig, prob_data_source = build_prob_figure(plasmy_fig)
	conservation_fig = plot_conservation_annotations(plasmy_fig, prob_data_source)

	# coverage filter
	# coverage_filter1, coverage_filter2= build_coverage_filter(plasmy_source)
	coverage_filter1 = build_coverage_filter(plasmy_source)

	# Build widgets
	search_input = build_search_input(plasmy_fig, label_source, line_source)
	clear_button = build_clear_button(label_source, line_source)

	# Layout figures and widgets
	layout_plots(
		plasmy_fig, 
		conservation_fig, 
		annotation_fig, 
		prob_fig, 
		coverage_filter1, 
		search_input, 
		clear_button
	)

#------------------------------------------------------------------------------
def acgt_color(base):
	# color = dict(A='#1f77b4', C='#9467bd', G='#2ca02c', T='#d62728')
	# color = dict(A='red', C='green', G='blue', T='black')
	color = dict(A=Set1[7][0], C=Set1[7][1], G=Set1[7][2], T=Set1[7][3], D=Set1[7][4], I=Set1[7][6])
	return color[base]

def plasmy_color(row):
	if row['A']>row['C'] and row['A']>row['G'] and row['A']>row['T'] and row['A']>row['D'] and row['A']>row['I']:
		return acgt_color('A')
	if row['C']>row['A'] and row['C']>row['G'] and row['C']>row['T'] and row['C']>row['D'] and row['C']>row['I']:
		return acgt_color('C')
	if row['G']>row['A'] and row['G']>row['C'] and row['G']>row['T'] and row['G']>row['D'] and row['G']>row['I']:
		return acgt_color('G')
	if row['T']>row['A'] and row['T']>row['C'] and row['T']>row['G'] and row['T']>row['D'] and row['T']>row['I']:
		return acgt_color('T')
	if row['D']>row['A'] and row['D']>row['C'] and row['D']>row['G'] and row['D']>row['T'] and row['D']>row['I']:
		return acgt_color('D')
	return acgt_color('I')

#------------------------------------------------------------------------------
def certainty(p):
	return 2.0 - sum([ -q*math.log2(q) for q in p if q>0] )

def plasmy_alpha(row):
	certainty_int = [certainty([0,0,0.5,0.5]), certainty([0,0,0.05,0.95])]   
	alpha_int = [0.4,1]
	min_alpha = 0.1
	h = certainty([row['A'],row['C'],row['G'],row['T'],row['D'],row['I']])
	return numpy.interp(h, certainty_int, alpha_int, left=min_alpha, right=1)

#------------------------------------------------------------------------------
# LAYOUT FIGURES AND WIDGETS
#------------------------------------------------------------------------------
def layout_plots(plasmy_fig, conservation_fig, annotation_fig, prob_fig, coverage_filter1, search_input, clear_button):
	acgt = figure(
		plot_width = DIM[0,1][0],
		plot_height = DIM[0,1][1],
		x_range = (0,5),
		y_range = (0.3,3),
		toolbar_location=None,
	)
	acgt.xgrid.grid_line_color = None
	acgt.ygrid.grid_line_color = None
	acgt.xaxis.visible = False
	acgt.xgrid.grid_line_color = None
	acgt.yaxis.visible = False
	acgt.ygrid.grid_line_color = None
	# acgt.min_border = 0
	acgt.outline_line_width = 1
	acgt.outline_line_alpha = 0.5
	acgt.outline_line_color = 'gray'

	source_A = ColumnDataSource(data=dict(
		x=[1,2,3,4,5,6],
		y=[1,1,1,1,1,1],
		text=['A','C','G','T','D','I'],
		text_color=[acgt_color('A'), acgt_color('C'), acgt_color('G'), acgt_color('T'), acgt_color('D'), acgt_color('I')],
	))
	lab_A = LabelSet(
		x='x',y='y',text='text',text_color='text_color',text_align='center',
		text_font_style = 'bold',
		source=source_A, level='glyph', render_mode='canvas')
	acgt.add_layout(lab_A)

	layout = column(
		row(
			column(plasmy_fig, conservation_fig, annotation_fig),
			column(prob_fig, acgt),
			column(widgetbox(coverage_filter1,clear_button,search_input, width=200)),
		),
	)
	print('Saved to', ARGS.output)
	output_file(ARGS.output, mode='inline', title='Heteroplasmy in %s' % ARGS.genome_name)
	# show(layout, browser="firefox")
	show(layout)

#------------------------------------------------------------------------------
# PLOT HETEROPLASMIC SITES
#------------------------------------------------------------------------------
def plot_heteroplasmies():
	# def assign_color(dist_to_neighbor):
	# 	for interval, color in NN_COLOR_SCHEME.items():
	# 		if interval[0] <= dist_to_neighbor <= interval[1]:
	# 			return color
	# 	return PALETTE_PLASMY[0]
	# PALETTE_PLASMY = Blues9[::-1]
	# NN_COLOR_SCHEME = {
	# 	(0,10) 			: PALETTE_PLASMY[8],
	# 	(11, 100) 		: PALETTE_PLASMY[7],
	# 	(101, 200) 		: PALETTE_PLASMY[6],
	# 	(201, np.inf)	: PALETTE_PLASMY[4]
	# }

	#---------------------------------------------------------------------------
	# Get data, build data source
	#---------------------------------------------------------------------------
	global MAX_X, VISIBLE_SAMPLE_RANGE

	plasmy_df = pandas.read_csv(ARGS.heteroplasmies)
	# plasmy_df['color'] = [ assign_color(value) for value in plasmy_df['d'] ]
	# plasmy_df['alpha'] = [ 1 for value in plasmy_df['total'] ]
	plasmy_df['color'] = [ plasmy_color(r[1]) for r in plasmy_df.iterrows() ]
	plasmy_df['alpha'] = [ plasmy_alpha(r[1]) for r in plasmy_df.iterrows() ]
	plasmy_df['alpha_original'] = [ plasmy_alpha(r[1]) for r in plasmy_df.iterrows() ]
		
	plasmy_source = ColumnDataSource(data = plasmy_df)
	if plasmy_df.max()['Coordinate'] > MAX_X:
		MAX_X = plasmy_df.max()['Coordinate']

	if VISIBLE_SAMPLE_RANGE[1] > plasmy_df['Sample'].max() + 1:
		VISIBLE_SAMPLE_RANGE = (VISIBLE_SAMPLE_RANGE[0], plasmy_df['Sample'].max() + 1)

	#---------------------------------------------------------------------------
	# Do the plotting
	#---------------------------------------------------------------------------
	p_hover = HoverTool(
		tooltips = [
			# ('Sample', '@Type, @Name'),
			('Sample', '@Name'),
			('Coordinate', '@Coordinate'),
			('Gene Product', '@GP'),
			('A', '@A{1.1111}'),
			('C', '@C{1.1111}'),
			('G', '@G{1.1111}'),
			('T', '@T{1.1111}'),
			('D', '@D{1.1111}'),
			('I', '@I{1.1111}'),
			('Coverage', '@total'),
			('NN distance', '@d'),
		],
		names = [ 'plasmy' ],
	)
	fig = figure(
		title='Heteroplasmy in %s' % ARGS.genome_name,
		plot_width = DIM[1,0][0],
		plot_height = DIM[1,0][1],
		tools=["xpan,ypan,xwheel_zoom,ywheel_zoom,box_zoom,undo,reset,save",p_hover],
		active_scroll="xwheel_zoom",
		y_range = VISIBLE_SAMPLE_RANGE,
		output_backend="webgl",
		logo=None,
		toolbar_location="above",
		# toolbar_sticky=False,
	)
	fig.xgrid.grid_line_color = None
	fig.xaxis.visible = False
	fig.ygrid.grid_line_color = None
	# fig.yaxis.visible = False

	person_id = plasmy_df['Sample']
	person_name = plasmy_df['Name']
	y_ticks_labels = { person_id[i] : person_name[i] for i in range(len(person_id)) }
	fig.axis.ticker = FixedTicker(ticks=person_id)
	fig.yaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    return labels[tick];
""" % y_ticks_labels )

	# plasmy_df = plasmy_df[plasmy_df['alpha'] == 1]
	plasmy_source = ColumnDataSource(data = plasmy_df)

	# fig.min_border = 0
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.circle(
		x = 'Coordinate',
		y = 'Sample',
		color = 'color',
		alpha = 'alpha',
		size = 6,
		name = 'plasmy',
		source = plasmy_source,
	)

	#---------------------------------------------------------------------------
	# Update some global variables (used by other plots)
	#---------------------------------------------------------------------------
	global HETEROPLASMY_PROBABILITIES

	g = plasmy_df[['Coordinate','Sample','A','C','G','T','D','I']].groupby('Coordinate')
	for gid in g.groups:
		rows = g.get_group(gid).iterrows()
		HETEROPLASMY_PROBABILITIES[gid] = [
			[r[1]['Sample'],r[1]['A'],r[1]['C'],r[1]['G'],r[1]['T'],r[1]['D'],r[1]['I']] for r in rows
		]

	return fig, plasmy_source

#------------------------------------------------------------------------------
# PLOT GENOME ANNOTATIONS
#------------------------------------------------------------------------------
def plot_genome_annotations(main_fig):
	global GENE_INTERVAL, MAX_X
	color_scheme = Accent8
	GENEPRODUCT_COLOR_SCHEME = dict(
		gene = color_scheme[0],
		exon = color_scheme[0],
		CDS = color_scheme[3],
		rRNA = color_scheme[1],
		tRNA = color_scheme[2],
		repeat_region = color_scheme[7],
	)
	#---------------------------------------------------------------------------
	# Get all entries
	#---------------------------------------------------------------------------
	entries = []
	locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )
	fill_alpha = dict(rRNA=1,tRNA=1,exon=0.25,gene=0.9,CDS=0.9,repeat_region=0.3)
	with open(ARGS.genome_annotations) as file:
		reader = csv.DictReader(file)
		for row in reader:
			if row['Direction'] in ['forward','reverse'] and row['Type'] in ['rRNA', 'tRNA', 'exon', 'gene', 'CDS', 'repeat_region']:
				if locale.atoi(row['Maximum']) > MAX_X:
					MAX_X = locale.atoi(row['Maximum'])
				y_coord = 1 if row['Direction']=='forward' else 0
				height =	0.2 if row['Type']=='repeat_region' else 1
				entries.append((
					locale.atoi(row['Minimum']),
					locale.atoi(row['Maximum']),
					row['Name'].strip(),
					row['Type'].strip(),
					row['# Intervals'],
					row['Direction'],
					y_coord,
					GENEPRODUCT_COLOR_SCHEME[row['Type']],
					height,
					fill_alpha[row['Type']],
					2 if row['Type']=='exon' else 1,				# linewidth
				))
				if row['Type'] != 'exon':
					int_min, int_max = locale.atoi(row['Minimum']), locale.atoi(row['Maximum'])
					name = row['Name'].split(' ')[0].strip()
					if name not in GENE_INTERVAL:
						GENE_INTERVAL[name.lower()] = [name, int_min, int_max, int_min, int_max]
					else:
						GENE_INTERVAL[name.lower()].append(int_min)
						GENE_INTERVAL[name.lower()].append(int_max)
						if int_min < GENE_INTERVAL[name][1]:
							GENE_INTERVAL[name.lower()][1] = int_min
						if int_max > GENE_INTERVAL[name][2]:
							GENE_INTERVAL[name.lower()][2] = int_max
	entries = sorted(entries)

	#---------------------------------------------------------------------------
	# Build data source
	#---------------------------------------------------------------------------
	source = ColumnDataSource(data=dict(
		x = [ (a[0]+a[1])//2 for a in entries ],
		y = [ a[6] for a in entries ],
		color = [ a[7] for a in entries ],
		height = [ a[8] for a in entries ],
		width = [ a[1]-a[0]+1 for a in entries ],
		lab = [ a[2] for a in entries ],
		min = [ a[0] for a in entries ],
		max = [ a[1] for a in entries ],
		dir = [ a[5] for a in entries ],
		fill_alpha = [ a[9] for a in entries ],
		line_alpha = [ a[9] for a in entries ],
		line_width = [ a[10] for a in entries ],
	))
	a_hover = HoverTool(
		tooltips = [
			('Name', '@lab'),
			('Location', '(@min, @max)'),
			('Direction', '@dir')
		],
		names = [ 'gene_product' ],
	)
	#---------------------------------------------------------------------------
	# Do the plotting
	#---------------------------------------------------------------------------
	if MAX_X > 1000000:
		MAX_X = 1000000
	main_fig.x_range = Range1d(-0.05*MAX_X, MAX_X*1.05)

	fig = figure(
		plot_width = DIM[2,0][0],
		plot_height= DIM[2,0][1],
		x_range = main_fig.x_range,
		y_range = (-3.5,2),
		tools=['reset,tap,xwheel_zoom',a_hover],
		toolbar_location=None,
		active_scroll="xwheel_zoom",
		# output_backend="webgl",
		logo=None,
	)
	fig.xgrid.grid_line_color = None
	fig.xaxis.axis_label_text_font_style = "normal"
	fig.xaxis.axis_label_text_font_size = "14pt"
	fig.xaxis[0].formatter = NumeralTickFormatter(format="0")
	fig.ygrid.grid_line_color = None
	fig.yaxis.visible = False
	# fig.min_border = 0
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'

	fig.rect(
		x = 'x',
		y = 'y',
		width = 'width',
		height = 'height',
		color = 'color',
		alpha = 'alpha',
		fill_alpha = 'fill_alpha',
		line_alpha = 'line_alpha',
		line_width = 'line_width',
		nonselection_color = 'color',
		width_units = 'data',
		height_units = 'data',
		name = 'gene_product',
		source = source
	)
	return fig

#------------------------------------------------------------------------------
# PLOT LABELLED RESULTS OF SEARCH
#------------------------------------------------------------------------------
def plot_search_result_annotations(annotation_fig):
	#---------------------------------------------------------------------------
	# Build data sources
	#---------------------------------------------------------------------------
	label_source = ColumnDataSource(data=dict(x=[],y=[],text=[]))
	line_source = ColumnDataSource(data=dict(xs=[],ys=[]))

	#---------------------------------------------------------------------------
	# Do the plotting
	#---------------------------------------------------------------------------
	annotation_fig.multi_line(
		xs = 'xs', ys = 'ys', line_color = 'navy', line_dash = 'dotted',
		line_alpha = 0.5, line_width = 2, source = line_source,
	)
	labels = LabelSet(
		x = 'x', y = 'y', text = 'text', text_align = 'center',
		level = 'glyph', render_mode = 'canvas', source = label_source,
	)
	annotation_fig.add_layout(labels)
	return label_source, line_source

#------------------------------------------------------------------------------
# PLOT CONSERVATION ANNOTATIONS
#------------------------------------------------------------------------------
def plot_conservation_annotations(main_fig, targeted_source):
	df = pandas.read_csv(ARGS.conserved_scores)
	#---------------------------------------------------------------------------
	# Build data source
	#---------------------------------------------------------------------------
	source = ColumnDataSource(data=dict(
		y = [0] * len(df),
		Coordinate = df['Coordinate'],
		Score = df['Score'],
	))
	source.callback = CustomJS(
		args=dict(
			targeted_source=targeted_source,
		),
		code="""
		var inds = cb_obj.selected['1d'].indices;
		var selected_data = cb_obj.data;
		var targeted_data = targeted_source.data;
		var prob = %s;
		var coord, items, i, v;
		targeted_data['y'] = [];
		targeted_data['left'] = [];
		targeted_data['right'] = [];
		targeted_data['height'] = [];
		targeted_data['color'] = [];
		samples = [];
		for (j=0; j<inds.length; j++) {
			coord = selected_data['Coordinate'][inds[j]]
			items = prob[coord];
			for (i=0; i<items.length; i++) {
				v = items[i];
				if (v[0] in samples) {
					u = samples[v[0]];
					samples[v[0]] = [u[0]+1, u[1]+v[1], u[2]+v[2], u[3]+v[3], u[4]+v[4], u[5]+v[5], u[6]+v[6]];
				} else {
					samples[v[0]] = [1, v[1], v[2], v[3], v[4], v[5], v[6]];
				}
			}
		}
		for (var s in samples) {
			if (samples.hasOwnProperty(s)) {
				u = samples[s];
				v = [u[1]/u[0], u[2]/u[0], u[3]/u[0], u[4]/u[0], u[5]/u[0], u[6]/u[0]];
				y = parseInt(s);
				Array.prototype.push.apply(targeted_data['y'], [y,y,y,y,y,y]);
				Array.prototype.push.apply(targeted_data['left'], [0,v[0],v[0]+v[1],v[0]+v[1]+v[2],v[0]+v[1]+v[2]+v[3],v[0]+v[1]+v[2]+v[3]+v[4]]);
				Array.prototype.push.apply(targeted_data['right'], [v[0],v[0]+v[1],v[0]+v[1]+v[2],v[0]+v[1]+v[2]+v[3],v[0]+v[1]+v[2]+v[3]+v[4],1]);
				Array.prototype.push.apply(targeted_data['height'], [0.9,0.9,0.9,0.9,0.9,0.9]);
				Array.prototype.push.apply(targeted_data['color'], ['%s','%s','%s','%s','%s','%s']);
			}
		}

		targeted_source.change.emit();
	""" % (HETEROPLASMY_PROBABILITIES,acgt_color('A'),acgt_color('C'),acgt_color('G'),acgt_color('T'),acgt_color('D'),acgt_color('I')))

	c_hover = HoverTool(
		tooltips = [
			('Coordinate', '@Coordinate'),
			('Conservation', '@Score{0,0.0000}'),
		],
		names = [ 'conserved' ],
	)

	#---------------------------------------------------------------------------
	# Do the plotting
	#---------------------------------------------------------------------------
	fig = figure(
		plot_width = DIM[0,0][0],
		plot_height= DIM[0,0][1],
		x_range = main_fig.x_range,
		tools=['tap,box_select,xwheel_zoom',c_hover],
		toolbar_location=None,
		active_scroll='xwheel_zoom',
		active_tap='tap',
		active_drag='box_select',
		logo = None,
		# webgl=True,
		output_backend="webgl",
	)
	fig.xgrid.grid_line_color = None
	fig.ygrid.grid_line_color = None
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.xaxis.visible = False
	fig.yaxis.visible = False
	# fig.min_border = 0

	# Reverse the color order so darkest has the highest value
	PALETTE_CONSERVATION_SCORE = YlOrRd9[::-1][2:]
	c_mapper = LinearColorMapper(PALETTE_CONSERVATION_SCORE, low=0, high=1)
	fig.square(
		x = 'Coordinate',
		y = 'y',
		color = {'field':'Score','transform':c_mapper},
		alpha = 1,
		size = 6,
		name = 'conserved',
		source = source,
	)
	return fig

#------------------------------------------------------------------------------
# This figure provides annotation of heteroplasmy probabilites across samples.
#------------------------------------------------------------------------------
def build_prob_figure(main_fig):
	fig = figure(
		plot_width =  DIM[1,1][0],
		plot_height = DIM[1,1][1],
		x_range = (0,1),
		y_range = main_fig.y_range,
		tools=[],
		toolbar_location=None,
	)
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.xaxis.visible = False
	fig.xgrid.grid_line_color = None
	fig.yaxis.visible = False
	fig.ygrid.grid_line_color = None
	# fig.min_border = 0

	prob_source = ColumnDataSource(data=dict(y=[],left=[],right=[],height=[],color=[]))
	prob_label_source = ColumnDataSource(data=dict(x=[],y=[],text=[]))


	#---------------------------------------------------------------------------
	# Do the plotting
	#---------------------------------------------------------------------------
	fig.hbar(y='y',left='left',right='right',color='color',height='height',source=prob_source)
	return fig, prob_source

#------------------------------------------------------------------------------
# Search provides zooming into a gene
#------------------------------------------------------------------------------
def build_search_input(main_fig, label_source, line_source):
	text = AutocompleteInput(
		title = 'Zoom into a specific gene',
		value = '',
		placeholder = 'Gene symbol',
		completions = [ v[0] for k,v in GENE_INTERVAL.items() ],
	)
	text.callback = CustomJS(
		args = dict(
			x_range = main_fig.x_range,
			label_source = label_source,
			line_source = line_source,
		),
		code="""
		var gene_symbol = cb_obj.value.toLowerCase();
		var interval = %s;
		var y = %s;
		var data = label_source.data;
		var data_line = line_source.data;
		var start, end, i;
		if (gene_symbol.length > 0) {
			if (gene_symbol in interval) {
				name = interval[gene_symbol][0];
				x_range.start = interval[gene_symbol][1] * 0.99;
				x_range.end = interval[gene_symbol][2] * 1.01;
				for (i=3; i<interval[gene_symbol].length; i += 2){
					start = interval[gene_symbol][i];
					end = interval[gene_symbol][i+1];
					data['x'].push((start+end)*0.5);
					data['y'].push(y);
					data['text'].push(name);
					data_line['xs'].push([start, end]);
					data_line['ys'].push([y,y]);
				}
				label_source.change.emit();
				line_source.change.emit();
			} else {
				;
			}
		}
		""" % (GENE_INTERVAL, -3))
	return text

#------------------------------------------------------------------------------
# THIS CLEARS GENE NAMES LABELED BY SEARCH
#------------------------------------------------------------------------------
def build_clear_button(label_source, line_source):
	button = Button(label='Clear gene labels', width=200)
	button_callback = CustomJS(
		args = dict(
			label_source = label_source,
			line_source = line_source,
		),
		code="""
		var data = label_source.data;
		var data_line = line_source.data;
		data['x'] = [];
		data['y'] = [];
		data['text'] = [];
		label_source.change.emit();

		data_line['xs'] = [];
		data_line['ys'] = [];
		line_source.change.emit();
	""")
	button.js_on_event(events.ButtonClick,button_callback)
	return button

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# COVERAGE FILTER
#------------------------------------------------------------------------------
def build_coverage_filter(plasmy_source):

	# round up to the next hundred
	def roundup(x):
		return int(math.ceil(x / 100.0)) * 100

	max_coverage = plasmy_source.data['total'].max()
	
	def slider_callback(source=plasmy_source, window=None):
		data = source.data
		slider_value = cb_obj.value
		total = data['total']
		alpha = data['alpha']
		alpha_original = data['alpha_original']
		for i in range(len(total)):
			alpha[i] = 0
			if total[i] < slider_value:
				alpha[i] = 0
			else:
				alpha[i] = alpha_original[i]

		source.change.emit()


	slider1 = Slider(start=0, end=roundup(max_coverage), value=0, step=100, title="Coverage 0-max", width = 200, callback=CustomJS.from_py_func(slider_callback))
	# slider2 = Slider(start=1000, end=roundup(max_coverage), value=1000, step=100, title="Coverage 1000x - max", width = 200, callback=CustomJS.from_py_func(slider_callback))

	# return slider1, slider2
	return slider1

#------------------------------------------------------------------------------

main()
