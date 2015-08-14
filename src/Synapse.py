"""
@author: Brett Settle
@Department: UCI Neurobiology and Behavioral Science
@Lab: Parker Lab
@Date: August 6, 2015
"""
import os,sys,inspect
from BioDocks import *
from pyqtgraph.Qt import QtCore, QtGui
from scipy.spatial import Delaunay
from pyqtgraph.dockarea import *
from collections import OrderedDict
import pyqtgraph.console

app = QtGui.QApplication([])

win = QtGui.QMainWindow()
win.resize(1700, 900)
dockArea = DockArea()
win.setCentralWidget(dockArea)
win.setWindowTitle('Main Window')
plotWidget = PlotWidget(viewBox=ROIViewBox(creatingROI=True, roiSnap=False))

units = {'Pixels': 166, 'Nanometers': 1}
unit_prefixes = {'Pixels': 'px', 'Nanometers': 'nm'}
unit = 'Nanometers'

ignore = {'Z Rejected'}

yAxis = plotWidget.getPlotItem().getAxis('left')
xAxis = plotWidget.getPlotItem().getAxis('bottom')
yAxis.setLabel(text='Y', units=unit)
yAxis.enableAutoSIPrefix(False)
xAxis.enableAutoSIPrefix(False)
xAxis.setLabel(text='X', units=unit)

channel_colors = [(255, 0, 0), (0, 255, 0)]
legend = plotWidget.addLegend()

Channels = []

class ActivePoint(QtCore.QPointF):
	def __init__(self, data):
		super(ActivePoint, self).__init__(data['Xc'], data['Yc'])
		self.data = data

	def __getitem__(self, item):
		if item in self.data:
			return self.data[item]
		else:
			return self.__dict__[item]


class Channel(pg.ScatterPlotItem):
	def __init__(self, name, points, **args):
		self.__name__ = name
		self.point_list = points
		base = {'brush': (255, 255, 255), 'size': 3, 'symbol': 'o', 'pen': None}
		base.update(args)
		super(Channel, self).__init__(x=np.array([p.x() for p in self.point_list]), y=np.array([p.y() for p in self.point_list]), **base)
		self.make_menu()
		self.colorDialog=QtGui.QColorDialog()
		self.colorDialog.currentColorChanged.connect(self.setBrush)

	def getCount(self):
		return len(self.point_list)

	def getPoints(self):
		return np.transpose([[p.x() for p in self.point_list], [p.y() for p in self.point_list]])

	def getCenter(self):
		if self.getCount() == 0:
			raise Exception('Cannot get center, no points in %s channel' % self.__name__)
		return np.average(self.getPoints(), 0)

	def make_menu(self):
		self.menu = QtGui.QMenu(self.__name__)
		self.menu.addAction(QtGui.QAction('Hide Item', self, triggered=lambda f: self.setVisible(not f), checkable=True))
		self.menu.addAction(QtGui.QAction("Change C&olor", self.menu, triggered =lambda :  self.colorDialog.open()))
		self.menu.addAction(QtGui.QAction('&Remove', self, triggered=lambda : self.getViewBox().removeItem(self)))

	def remove(self):
		plotWidget.removeItem(self)
		Channels.remove(self)

def clear():
	global Channels
	for ch in Channels:
		legend.removeItem(ch.__name__)
	Channels = []
	plotWidget.clear()

def import_channels(filename=''):
	clear()
	if filename == '':
		filename = getFilename('Select the text file of the channels to import')
	if filename == '':
		return
	data = file_to_arr(filename)
	data['Channel Name'] = data['Channel Name'].astype(str)
	data['Xc'] /= units[unit]
	data['Yc'] /= units[unit]
	channel_names = {str(i) for i in data['Channel Name']} - ignore
	assert len(channel_names) == 2, 'Must provide only 2 channels, channels are %s' % channel_names
	pts = [ActivePoint({k: data[k][i] for k in data}) for i in range(len(data['Channel Name']))]

	for i, ch in enumerate(channel_names):
		item = Channel(name=ch, points=[p for p in pts if p['Channel Name'] == ch], brush=channel_colors[i])
		plotWidget.add(item, name=ch)
		Channels.append(item)
		legend.addItem(item, ch)

def displayData():
	synapseWidget.setData(sorted([roi.synapse_data for roi in plotWidget.items() if isinstance(roi, Freehand) and hasattr(roi, 'synapse_data')], key=lambda f: f['ROI #']))

def analyze_roi(roi):
	channels = []
	for ch in Channels:
		pts_in = []
		for syn_pt in ch.point_list:
			if roi.contains(syn_pt):
				pts_in.append(syn_pt)
		channels.append(Channel(ch.__name__, pts_in))

	roi.synapse_data = OrderedDict([('ROI #', roi.id), ('Mean Distance (%s)' % unit_prefixes[unit], 0), ('%s N' % Channels[0].__name__, 0), \
	('%s N' % Channels[1].__name__, 0), ('%s Area (%s^2)' % (Channels[0].__name__, unit_prefixes[unit]), 0), ('%s Area (%s^2)' % (Channels[1].__name__, unit_prefixes[unit]), 0)])

	for i, ch in enumerate(channels):
		roi.synapse_data['%s N' % ch.__name__] = ch.getCount()
		if ch.getCount() >= 3:
			roi.synapse_data['%s Area (%s^2)' % (ch.__name__, unit_prefixes[unit])] = concaveArea(ch.getPoints(), channel_colors[i])
		else:
			print('Cannot get area of %s in roi %d with %d points' % (ch.__name__, roi.id, ch.getCount()))

	if hasattr(roi, 'mean_line'):
		plotWidget.removeItem(roi.mean_line)

	if all([ch.getCount() > 0 for ch in channels]):
		roi.synapse_data['Mean Distance (%s)' % unit_prefixes[unit]] = np.linalg.norm(channels[1].getCenter() - channels[0].getCenter())
		roi.mean_line = pg.PlotDataItem(np.array([channels[1].getCenter(), channels[0].getCenter()]), symbol='d')
		roi.mean_line.setParentItem(roi)
		roi.mean_line.setVisible(False)
	else:
		del roi.synapse_data
		print('Must select exactly 2 channels to calculate distance. Ignoring ROI %d' % roi.id)

	displayData()

def distance(ptA, ptB):
	return np.linalg.norm(np.subtract(ptA, ptB))

def order_walls(walls):
	new_wall = walls.pop(0)
	while walls:
		add = [wall for wall in walls if new_wall[-1] in wall][0]
		walls.remove(add)
		add.remove(new_wall[-1])
		new_wall.extend(add)
	return new_wall

def getTriangleArea(A, B, C):
	return .5 * abs(A[0]*(B[1] - C[1]) + B[0]*(C[1] - A[1]) + C[0]*(A[1] - B[1]))

def concaveArea(points, color):
	tri = Delaunay(points)
	outerwalls = tri.convex_hull.tolist()
	outerwalls = order_walls(outerwalls)
	verts = tri.vertices.tolist()
	change = False
	i = 0
	while i < len(outerwalls) - 1:
		at = outerwalls[i]
		next = outerwalls[i + 1]
		outer_dist = distance(points[at], points[next])
		inner = None
		for t in verts:
			inners = set(t) ^ {at, next}
			if len(inners) == 1 and len(set(outerwalls) & set(t)) == 2:
				inner = inners.pop()
				break
		if inner != None and outer_dist > distance(points[at], points[inner]):
			outerwalls.insert(i+1, inner)
			change = True
			verts.remove(t)
			i += 1
		i += 1
		if i >= len(outerwalls) - 1 and change:
			change = False
			i = 0
	pts = np.array([points[i] for i in outerwalls])
	#for vs in verts:
	#	plotWidget.addItem(pg.PlotDataItem(np.array([points[i] for i in (vs + [vs[0]])]), pen=color))
	return sum(map(lambda vs: getTriangleArea(*[points[i] for i in vs]), verts))

def show_line(roi, hover):
	if hasattr(roi, 'mean_line'):
		roi.mean_line.setVisible(hover)

def connect_roi(roi):
	roi.sigChanged.connect(analyze_roi)
	roi.sigRemoved.connect(lambda : displayData())
	roi.sigHoverChanged.connect(show_line)
	analyze_roi(roi)

menu = win.menuBar()

fileMenu = QtGui.QMenu('&File', menu)
fileMenu.addAction(QtGui.QAction('&Import Channels', fileMenu, triggered = lambda : import_channels()))
fileMenu.addAction(QtGui.QAction('&Close', fileMenu, triggered = win.close))
menu.addMenu(fileMenu)
plotWidget.load_file = import_channels
plotDock = WidgetDock(name='Plot Dock', size=(500, 400), widget=plotWidget)
dockArea.addDock(plotDock)

synapseWidget = DataWidget()
synapseWidget.setFormat("%3.3f")
synapseDock = Dock(name='MeanXY Distances', widget=synapseWidget)
dockArea.addDock(synapseDock, 'right', plotDock)

plotWidget.getViewBox().roiCreated.connect(connect_roi)

win.show()
sys.exit(app.exec_())
