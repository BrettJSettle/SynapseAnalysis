"""
@author: Brett Settle
@Department: UCI Neurobiology and Behavioral Science
@Lab: Parker Lab
@Date: August 6, 2015
"""
import os,sys
from BioDocks import *
from pyqtgraph.dockarea import *
from scipy.spatial import ConvexHull
from collections import OrderedDict
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from Channels import *
from ClusterMath import *

app = QApplication([])

win = DockWindow(addMenu=False)
win.resize(1700, 900)
win.setWindowTitle('Main Window')
plotWidget = PlotWidget(viewBox=ROIViewBox(creatingROI=True, roiSnap=False))

# data is loaded in nanometers, divided by # according to units
units = {'Pixels': 166, 'Nanometers': 1}
unit_prefixes = {'Pixels': 'px', 'Nanometers': 'nm'}
unit = 'Nanometers'

yAxis = plotWidget.getPlotItem().getAxis('left')
xAxis = plotWidget.getPlotItem().getAxis('bottom')
yAxis.setLabel(text='Y', units=unit)
yAxis.enableAutoSIPrefix(False)
xAxis.enableAutoSIPrefix(False)
xAxis.setLabel(text='X', units=unit)

legend = plotWidget.addLegend()

colors = ((255, 0, 0), (0, 255, 0))
color_dict = {'atto488': colors[0], 'Alexa647': colors[1]}
ignore = {"Z Rejected"}

Channels = []
empty_channel = Channel('Empty', [], (1, 1, 1))

def displayData():
	dataWidget.setData(sorted([roi.synapse_data for roi in plotWidget.items() if \
		isinstance(roi, Freehand) and hasattr(roi, 'synapse_data')], key=lambda f: f['ROI #']))
	dataWidget.changeFormat('%.3f')

def subchannels_in_roi(roi):
	channels = []
	for ch in Channels:
		pts_in = []
		for syn_pt in ch.pts:
			if roi.contains(QPointF(syn_pt[0], syn_pt[1])):
				pts_in.append(syn_pt)
		channels.append(Channel(ch.__name__, pts_in, ch.color()))
	return channels

def analyze_roi(roi):
	channels = subchannels_in_roi(roi)
	roi.synapse_data = OrderedDict([('ROI #', roi.id), ('Mean Distance (%s)' % unit_prefixes[unit], 0), ('%s N' % Channels[0].__name__, 0), \
	('%s N' % Channels[1].__name__, 0), ('%s Volume (%s^3)' % (Channels[0].__name__, unit_prefixes[unit]), 0), ('%s Volume (%s^3)' % (Channels[1].__name__, unit_prefixes[unit]), 0)])

	for i, ch in enumerate(channels):
		roi.synapse_data['%s N' % ch.__name__] = ch.getCount()
		if ch.getCount() >= 4:
			roi.synapse_data['%s Volume (%s^3)' % (ch.__name__, unit_prefixes[unit])] = convex_volume(ch.getPoints(True))
		else:
			print('Cannot get Volume of %s in roi %d with %d points' % (ch.__name__, roi.id, ch.getCount()))

	if hasattr(roi, 'mean_line'):
		plotWidget.removeItem(roi.mean_line)

	if all([ch.getCount() > 0 for ch in channels]):
		roi.synapse_data['Mean Distance (%s)' % unit_prefixes[unit]] = np.linalg.norm(channels[1].getCenter(z=True) - channels[0].getCenter(z=True))
		roi.mean_line = pg.PlotDataItem(np.array([channels[1].getCenter(), channels[0].getCenter()]), symbol='d')
		roi.mean_line.setParentItem(roi)
		roi.mean_line.setVisible(False)
	else:
		del roi.synapse_data
		print('Must select exactly 2 channels to calculate distance. Ignoring ROI %d' % roi.id)
	displayData()

def plotROIChannels(roi):
	if not synapseDock.isVisible():
		synapseDock.float()
		synapseDock.window().setGeometry(20, 100, 500, 500)
	ch1, ch2 = subchannels_in_roi(roi)
	if hasattr(synapseWidget, 'synapse'):
		synapseWidget.synapse.setChannels(ch1, ch2)
	else:
		synapseWidget.synapse = Synapse(ch1, ch2)
		synapseWidget.addItem(synapseWidget.synapse, name=roi.__name__)
	cen = np.average(synapseWidget.synapse.pos, axis=0)
	d = max([np.linalg.norm(np.subtract(p, cen)) for p in synapseWidget.synapse.pos])
	synapseWidget.moveTo(cen, distance = 2 * d)

def roiCreated(roi):
	roi.sigChanged.connect(analyze_roi)
	roi.sigRemoved.connect(lambda : displayData())
	roi.sigHoverChanged.connect(lambda r, h: r.mean_line.setVisible(h) if hasattr(r, 'mean_line') else None)
	roi.sigClicked.connect(plotROIChannels)
	analyze_roi(roi)

def open_file(filename=''):
	if filename == '':
		filename = getFilename(filter='Text Files (*.txt)')
	clear()
	data = importFile(filename)
	data = {d[0]: d[1:] for d in np.transpose(data)}
	for k in data:
		if k != 'Channel Name':
			data[k] = data[k].astype(float)
	print('Gathering channels...')
	names = set(data['Channel Name'].astype(str)) - ignore
	print('Channels Found: %s' % ', '.join(names))

	data['Xc'] /= units[unit]
	data['Yc'] /= units[unit]
	data['Zc'] /= units[unit]

	global Channels
	Channels = []
	plotWidget.clear()
	pts = [ActivePoint(data={k: data[k][i] for k in data}) for i in range(len(data['Channel Name']))]
	for i, n in enumerate(names):
		if n in color_dict:
			color = color_dict[n]
		else:
			color = colors[i]
		Channels.append(Channel(n, [p for p in pts if p['Channel Name'] == n], color))
		plotWidget.addItem(Channels[-1])
		legend.addItem(Channels[-1], n)
	show_ch1.setText(Channels[0].__name__)
	show_ch2.setText(Channels[1].__name__)
	ch1_mesh.setText(Channels[0].__name__)
	ch2_mesh.setText(Channels[1].__name__)

def show_mesh(i):
	if i == None:
		synapseWidget.synapse.mesh.setVisible(False)
		synapseWidget.synapse.update()
		return
	else:
		ps = synapseWidget.synapse.channels[Channels[i].__name__].getPoints(True)
		c = synapseWidget.synapse.channels[Channels[i].__name__].color()
	if len(ps) == 0:
		return
	tri = ConvexHull(ps)
	synapseWidget.synapse.mesh.setColor(c)
	synapseWidget.synapse.mesh.setMeshData(vertexes=ps, faces=tri.simplices)
	synapseWidget.synapse.mesh.meshDataChanged()
	synapseWidget.synapse.mesh.setVisible(True)
	synapseWidget.synapse.update()

def hide_channel(i):
	global Channels
	Channels[0].setVisible(True)
	Channels[1].setVisible(True)
	if i != None:
		Channels[i].setVisible(False)

def clear():
	global Channels
	for ch in Channels:
		legend.removeItem(ch.__name__)
	Channels = []
	plotWidget.clear()
	for i in plotWidget.getViewBox().addedItems:
		if isinstance(i, Freehand):
			i.delete()

menu = win.menuBar()

fileMenu = menu.addMenu('&File')
fileMenu.addAction(QAction('&Import Channels', fileMenu, triggered = lambda : open_file()))
fileMenu.addAction(QAction('&Close', fileMenu, triggered = win.close))

plotWidget.getViewBox().roiCreated.connect(roiCreated)
plotWidget.load_file = open_file

opsFrame = QWidget()
layout = QGridLayout(opsFrame)
show_all = QRadioButton('Show all', checked=True)
show_all.pressed.connect(lambda : hide_channel(None))
show_ch1 = QRadioButton('Channel 1')
show_ch1.pressed.connect(lambda : hide_channel(1))
show_ch2 = QRadioButton('Channel 2')
show_ch2.pressed.connect(lambda : hide_channel(0))
layout.addWidget(show_all, 0, 0)
layout.addWidget(show_ch1, 0, 1)
layout.addWidget(show_ch2, 0, 2)
layout.addItem(QSpacerItem(400, 20), 0, 3, 1, 8)

plotWidget.__name__ = '2D Plotted Channels'
plotDock = win.addWidget(plotWidget, size=(400, 500))
plotDock.addWidget(opsFrame)

synapseFrame = QWidget()
layout = QGridLayout(synapseFrame)
synapseWidget = Plot3DWidget()
synapseWidget.load_file = open_file
layout.addWidget(synapseWidget, 0, 0, 6, 6)
export_synapse = QPushButton('Export Coordinates')
export_synapse.pressed.connect(lambda : export_arr(synapseWidget.synapse.pos, header='X\tY\tZ'))
layout.addWidget(export_synapse, 6, 0)
no_mesh = QRadioButton('No Mesh', checked=True)
no_mesh.pressed.connect(lambda : show_mesh(None))
ch1_mesh = QRadioButton('Channel 1')
ch1_mesh.pressed.connect(lambda : show_mesh(0))
ch2_mesh = QRadioButton('Channel 2')
ch2_mesh.pressed.connect(lambda : show_mesh(1))
layout.addWidget(no_mesh, 6, 1)
layout.addWidget(ch1_mesh, 6, 2)
layout.addWidget(ch2_mesh, 6, 3)
layout.addWidget(QLabel(\
'''
ROI Plot 3D Widget Controls
	Arrow Keys or Left click and drag to rotate camera
	Middle click and drag to pan
	Scroll mouse wheel to zoom
	Right Click for plotted item options
'''), 7, 0, 2, 3)

synapseDock = win.addWidget(size=(300, 100), widget=synapseFrame)
synapseDock.hide()
plotDock.window().setGeometry(340, 100, 1400, 800)

dataWidget = DataWidget()
win.addWidget(dataWidget, where=('right', plotDock), size=(100, 500))
win.closeEvent = lambda f: app.exit()

win.show()
sys.exit(app.exec_())
